import sys
import itertools
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import string
import json
import collections
import os
import xml.etree.ElementTree as ET
import copy
import math

revcomp = None
if sys.version_info[0] >= 3:
    revcomp = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
else:
    import itertools.izip as zip
    revcomp = string.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')

VALID_BASES = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'n', 'N'}

I5 = 'i5'
I7 = 'i7'
R1 = 'r1'
R2 = 'r2'
BC_READ = 'read'
BC_START = 'start'
BC_END = 'end'
BC_WHITELIST = 'whitelist'
BC_CORRECTION_MAP = 'correction_map'

_accepted_read_keys = {I5, I7, R1, R2}
_accepted_barcode_properties = {BC_READ, BC_START, BC_END, BC_WHITELIST}
_required_barcode_properties = {BC_READ, BC_START, BC_END}
_reserved_keys = {'r1_name', 'r2_name', 'r1_seq', 'r2_seq', 'r1_qual', 'r2_qual'}


def get_run_info(flow_cell_path):
    """
    Helper function to get some info about the sequencing runs.

    Args:
        flow_cell_path (str): Path to BCL directory for run.

    Returns:
        dict: basic statistics about run, like date, instrument, number of lanes, flowcell ID, read lengths, etc.
    """
    run_stats = {}

    bcl_run_info = os.path.join(flow_cell_path, 'RunInfo.xml')

    if not os.path.exists(bcl_run_info):
        raise ValueError('BCL RunInfo.xml not found for specified flowcell: %s' % bcl_run_info)

    # Get the flowcell info
    flowcells = list(ET.parse(bcl_run_info).getroot().iter('Flowcell'))
    if len(flowcells) > 1:
        raise ValueError('BCL RunInfo contains %s flowcell declarations. Expected 1.' % len(flowcells))

    flowcell = flowcells[0].text.split('-')[-1]
    run_stats['flow_cell_id'] = flowcell

    date = list(ET.parse(bcl_run_info).getroot().iter('Date'))[0].text
    run_stats['date'] = date

    instrument = list(ET.parse(bcl_run_info).getroot().iter('Instrument'))[0].text
    run_stats['instrument'] = instrument

    flowcell_lanes = int(list(ET.parse(bcl_run_info).getroot().iter('FlowcellLayout'))[0].attrib['LaneCount'])
    run_stats['lanes'] = flowcell_lanes

    # Get the index read lengths
    read_info = list(list(ET.parse(bcl_run_info).getroot().iter('Reads'))[0].iter())

    if not read_info:
        raise ValueError('Read info not found in BCL file: %s' % bcl_run_info)

    for read in read_info:
        attributes = read.attrib
        if not 'IsIndexedRead' in attributes:
            continue

        if attributes['IsIndexedRead'] == 'N' and 'r1_length' not in run_stats:
            run_stats['r1_length'] = int(attributes['NumCycles'])
            continue
        if attributes['IsIndexedRead'] == 'Y' and 'p7_index_length' not in run_stats:
            run_stats['p7_index_length'] = int(attributes['NumCycles'])
            continue
        if attributes['IsIndexedRead'] == 'N' and 'r2_length' not in run_stats:
            run_stats['r2_length'] = int(attributes['NumCycles'])
            continue
        if attributes['IsIndexedRead'] == 'Y' and 'p5_index_length' not in run_stats:
            run_stats['p5_index_length'] = int(attributes['NumCycles'])
            continue

    return run_stats


def correct_barcode(barcode, mismatch_map):
    """
    Correct an observed raw barcode to one of a list of whitelists of mismatches.
    Args:
            barcode (string): barcode sequence to be corrected
            mismatch_map (list of dict dict): list of dict of mismatched sequences to real sequences
    Returns:
            string: corrected barcodes or None if barcode not correctable.
    """
    for mismatch_whitelist in mismatch_map:
        corrected = mismatch_whitelist.get(barcode, None)

        if corrected:
            return corrected

    return None


def generate_mismatches(sequence, num_mismatches, allow_n=True):
    """
    Generate a list of mimatched sequences to a given sequence. Must only contain ATGC.
    This is heavily based on a biostars answer.
    Args:
        sequence (str): The sequence must contain only A, T, G, and C
        num_mismatches (int): number of mismatches to generate sequences for
        allow_n (bool): True to allow N bases and False if not
    
    Returns:
        list: list of all specified mismatches
    """
    letters = 'ACGT'

    if allow_n:
        letters += 'N'

    sequence = sequence.upper()
    mismatches = []

    for locs in itertools.combinations(range(len(sequence)), num_mismatches):
        sequence_list = [[char] for char in sequence]
        for loc in locs:
            orig_char = sequence[loc]
            sequence_list[loc] = [l for l in letters if l != orig_char]

        for poss in itertools.product(*sequence_list):
            mismatches.append(''.join(poss))

    return mismatches


def construct_mismatch_to_whitelist_map(whitelist, edit_distance, allow_n=True):
    """
    Constructs a precomputed set of all mimatches within a specified edit distance and the barcode whitelist.
    Args:
        whitelist (set of str): set of whitelist sequences
        edit_distance (int): max edit distance to consider
        allow_n (bool): True to allow N bases and False if not
    Returns:
        dict: mapping of mismatched sequences to their whitelist sequences
    """

    mismatch_to_whitelist_map = [None] * (edit_distance + 1)

    mismatch_to_whitelist_map[0] = {k: k for k in whitelist}

    conflicting_mismatches = []  # tracks conflicts where mismatches map to different sequences

    # Doesn't really matter as correction function will never see it,
    # but exclude any perfect matches to actual seqs by mismatches
    conflicting_mismatches.extend(list(whitelist))

    for mismatch_count in range(1, edit_distance + 1):
        mismatch_to_whitelist_map[mismatch_count] = {}

        for sequence in whitelist:
            sequence = sequence.upper()

            # Generate all possible mismatches in range
            mismatches = generate_mismatches(sequence, num_mismatches=mismatch_count, allow_n=allow_n)

            # Construct a mapping to the intended sequences
            for mismatch in mismatches:
                # Check for conflict with existing sequence and track if so
                if mismatch in mismatch_to_whitelist_map[mismatch_count]:
                    conflicting_mismatches.append(mismatch)
                mismatch_to_whitelist_map[mismatch_count][mismatch] = sequence

        # Go back and remove any conflicting mismatches
        for mismatch in set(conflicting_mismatches):
            if mismatch in mismatch_to_whitelist_map[mismatch_count]:
                del mismatch_to_whitelist_map[mismatch_count][mismatch]

    return mismatch_to_whitelist_map

def is_dna(seq):
    """
    Tests for valid DNA sequence w.r.t. the functionality of this tool (A, T, G, C, N).

    Args:
        seq (str): DNA sequence to be tested.

    Returns:
        bool: True if valid DNA and False otherwise.
    """
    if not isinstance(seq, str):
        raise ValueError('Argument must be a string, found %s.' % type(seq))

    for base in seq:
        if base not in VALID_BASES:
            return False
    return True

def hamming_distance(str1, str2, capdistance=None):
    """Count the # of differences between equal length strings str1 and str2. Max diff can be capped using capdistance to short-circuit full calculation when only care about a range of distances."""
    diffs = 0

    if len(str1) != len(str2):
        raise ValueError('str1 and str2 must be equal lengths, but found %s and %s.' % (len(str1), len(str2)))

    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
            if capdistance is not None and diffs >= capdistance:
                return diffs
    return diffs

def reverse_complement(seq):
    return seq.translate(revcomp)[::-1]

def open_file(f, mode='rt'):
    if f.endswith('.gz'):
        return gzip.open(f, mode)
    else:
        return open(f, mode)

def validate_barcode_spec(spec):
    """
    Given a spec, validate that it is correct.

    Args:
        dict: Spec object

    Returns:
        bool: (True, None) or (False, error)
    """

    if not isinstance(spec, dict):
        return (False, 'Spec must be a dictionary and %s found.' % type(spec))

    if len(spec) == 0:
        return (False, 'Spec must not be empty dict.')

    for k in spec:
        if not isinstance(spec[k], dict):
            return (False, '%s key entry in dictionary must be a dictionary specifying the properties for a barcode, but %s object found.' % (k, type(k)))

        if k in _reserved_keys:
            return (False, 'Entry %s is shares the name of a key that we already return by default. Please change this name to something else.' % k)

        properties_seen = set()
        barcode_start = None
        barcode_end = None
        for property_key in spec[k]:
            properties_seen.add(property_key)
            bc_property = spec[k][property_key]

            if property_key == BC_START:
                barcode_start = spec[k][BC_START]
            elif property_key == BC_END:
                barcode_end = spec[k][BC_END]

            if (property_key == BC_START or property_key == BC_END) and not isinstance(bc_property, int):
                return (False, '%s property in entry for %s must be an int and %s found.' % (property_key, k, type(bc_property)))
            
            elif property_key == BC_READ and (not isinstance(property_key, str) or bc_property not in _accepted_read_keys):
                return (False, '%s property in entry for %s must be a string and %s found.' % (property_key, k, type(bc_property)))

            elif property_key == BC_WHITELIST:
                if isinstance(bc_property, str):
                    if not os.path.exists(bc_property):
                        return (False, '%s property in entry for %s is supposed to be a string or a set object. Found string, but no file is found with the name %s.' % (property_key, k, bc_property))
                elif not isinstance(bc_property, set) and not isinstance(bc_property, list):
                    return (False, '%s property in entry for %s is supposed to be a string or a set/list object.' % (property_key, k))
                else:
                    # Is set or list, check elements
                    if len(bc_property) == 0:
                        return (False, 'Set or list found for whitelist provided by %s, %s, but is empty.' % (property_key, bc_property))
                    for entry in bc_property:
                        if not isinstance(entry, str):
                            return (False, 'Entry in whitelist for %s, %s, is not a string.' % (k, entry))
                        if not is_dna(entry):
                            return (False, 'Entry in whitelist for %s, %s, is not a valid DNA sequence.' % (k, entry))
                        
            elif property_key not in _accepted_barcode_properties:
                return (False, 'Invalid entry found in spec: %s. See documentation for allowed keys.' % (property_key))
        
        if barcode_end < barcode_start:
            return (False, 'End specified for %s is less than the start (start: %s, end: %s).' % (k, barcode_start, barcode_end))

        for prop in _required_barcode_properties:
            if prop not in properties_seen:
                return (False, 'Required property "%s" not seen in spec for entry %s' % (prop, k))
    return True, None

def load_barcode_spec(f):
    """
    Loads barcodes spec from JSON file.

    Args:
        str or file: name of file or a file handle

    Returns:
        dict: object stored in JSON file.
    """
    if hasattr(f, 'read'):
        file_handle = f
    else:
        file_handle = open_file(f)
    
    try:
        spec = json.load(file_handle)
    except json.decoder.JSONDecodeError as e:
        raise ValueError('Spec JSON file could not be loaded from file, please see documentation for formatting details: %s' % e)
    
    valid, error = validate_barcode_spec(spec)
    if not valid:
        raise ValueError(error)
    
    return spec


def valid_whitelist(whitelist):
    observed_lengths = set()
    for seq in whitelist:
        if not is_dna(seq):
            return (False, 'Whitelist entry %s is not a DNA sequence and only DNA sequences are allowed (ATGCN).' % seq)
        observed_lengths.add(len(seq))
    
    if len(observed_lengths) > 1:
        return (False, 'Whitelist has barcodes of variable lengths %s' % (', '.join(list(observed_lengths))))

    return True, None

def load_whitelist(whitelist):
    if not os.path.exists(whitelist):
        raise ValueError('Specified whitelist file does not exist %s' % whitelist)

    whitelist = [line.strip() for line in open_file(whitelist)]
    whitelist = set(whitelist)
    valid, error = valid_whitelist(whitelist)
    if valid:
        return whitelist
    else:
        raise ValueError(error)

class BarcodeCorrecter:
    def __init__(self, whitelist, edit_distance=1):
        if not isinstance(whitelist, set) and not isinstance(whitelist, list):
            raise ValueError('Whitelist argument must be a list or a set.')
        
        if isinstance(whitelist, str):
            # this already validates so no need to do again
            whitelist = load_whitelist(whitelist)
        else:
            if isinstance(whitelist, list):
                whitelist = set(whitelist)
            
            valid, error = valid_whitelist(whitelist)
            if not valid:
                raise ValueError(error)
      
        if not isinstance(edit_distance, int) and edit_distance > 0:
            raise ValueError('Edit distance must be a positive integer and got %s.' % edit_distance)

        self.whitelist = whitelist
        self.edit_distance = edit_distance
        self.mismatch_map = construct_mismatch_to_whitelist_map(self.whitelist, self.edit_distance)

    def correct(self, seq):
        return correct_barcode(seq, self.mismatch_map)

    def get_min_hamming(self, n=1):
        """
        Returns a list of the minimum N hamming distances observed.
        """
        if not isinstance(n, int):
            raise ValueError('n argument must be an int, but %s found.' % n)
        if n <= 0:
             raise ValueError('n argument must be a positive int, but %s found.' % n)

        hamming_distances = []
        observed_pairs = set()
        for barcode in self.whitelist:
            for other_barcode in self.whitelist:
                pair = tuple(sorted([barcode, other_barcode]))

                if barcode == other_barcode or pair in observed_pairs:
                    continue
                observed_pairs.add(pair)
                hamming_distances.append(hamming_distance(barcode, other_barcode))

        return sorted(hamming_distances)[0:n]

    def get_barcode_length(self):
        """
        Returns the length of the barcodes in whitelist.
        """
        return len(list(self.whitelist)[0])

def _get_index_coords(r1_name):
    """
    Helper functions to get index read start and end coords given R1 read name.
    """
    i7_start_coord = r1_name.rfind(':') + 1

    # Handle not having any valid index seqs in read name
    if not is_dna(r1_name[i7_start_coord:].replace('+', '')):
        return None, None, None, None

    indices = r1_name[i7_start_coord:]

    if '+' in indices:
        split = r1_name.rfind('+')
        i5_start_coord = split + 1
        i5_end_coord = len(r1_name)
        i7_end_coord = split
    else:
        i5_start_coord = None
        i5_end_coord = None
        i7_end_coord = len(r1_name)
    
    return i7_start_coord, i7_end_coord, i5_start_coord, i5_end_coord


def _validate_barcode_read_pair(read_seq, bc_end):
    if bc_end > len(read_seq):
        raise ValueError('Requested region extends beyond end of %s bp long read.' % (len(read_seq)))

def parse_fastq_barcodes(r1, r2=None, spec=None, reverse_i5=False, edit_distance=2):
    spec = copy.deepcopy(spec)

    if not spec:
        raise ValueError('spec is a required argument and may not be None.')
    valid_spec, error = validate_barcode_spec(spec)

    if not valid_spec:
        raise ValueError(error)

    # Validate requested reads    
    requested_reads = set([spec[barcode][BC_READ] for barcode in spec])
    if not r2 and R2 in requested_reads:
        raise ValueError('%s requested but r2 not provided as an argument to parser...' % R2)

    # Spec is valid and r1/r2 valid, set up any whitelists
    for barcode in spec:
        if BC_WHITELIST in spec[barcode]:
            spec[barcode][BC_CORRECTION_MAP] = BarcodeCorrecter(spec[barcode][BC_WHITELIST], edit_distance=edit_distance)

            # Also, make sure requested length is in-line with their specified start and end 
            barcode_length = (spec[barcode][BC_END] - spec[barcode][BC_START] + 1)
            whitelist_length = spec[barcode][BC_CORRECTION_MAP].get_barcode_length()
            if barcode_length != whitelist_length:
                raise ValueError('The barcodes specified in your whitelist are not the same length as the requested region for barcode %s, %s to %s (length: %s; whitelist length: %s).' % (barcode, BC_START, BC_END, barcode_length, whitelist_length))

    # Now set up file handles
    if r1:
        if not hasattr(r1, 'read'):
            r1_handle = FastqGeneralIterator(open_file(r1))
        else:
            r1_handle = FastqGeneralIterator(r1)
    if r2:
        if not hasattr(r1, 'read'):
            r2_handle = FastqGeneralIterator(open_file(r2))
        else:
            r2_handle = FastqGeneralIterator(r2)

    if r1 is not None and r2 is not None:
        fastq_iterator = zip(r1_handle, r2_handle)
        paired_end = True
    else:
        fastq_iterator = r1_handle
        paired_end = False

    # Iterate
    barcodes_dict = dict()

    i5_start = None
    i5_end = None
    i7_start = None
    i7_end = None
    read_lengths_checked = False

    for r in fastq_iterator:
        # Get seqs and store
        if paired_end:
            (r1_name, r1_seq, r1_qual), (r2_name, r2_seq, r2_qual) = r
            barcodes_dict['r1_name'] = r1_name
            barcodes_dict['r1_seq'] = r1_seq
            barcodes_dict['r1_qual'] = r1_qual

            barcodes_dict['r2_name'] = r2_name
            barcodes_dict['r2_seq'] = r2_seq
            barcodes_dict['r2_qual'] = r2_qual            
        else:
            r1_name, r1_seq, r1_qual = r
            barcodes_dict['r1_name'] = r1_name
            barcodes_dict['r1_seq'] = r1_seq
            barcodes_dict['r1_qual'] = r1_qual

        # Find index coords
        if (i7_start is None and i7_end is None) and (I7 in requested_reads or I5 in requested_reads):
            i7_start, i7_end, i5_start, i5_end = _get_index_coords(r1_name)
            if (i5_start is None or i5_end is None) and I5 in requested_reads:
                raise ValueError('Spec requests I5 read in barcode %s, but i5 index not found in r1 name %s.' % (barcode, r1_name))
            if (i7_start is None or i7_end is None) and I7 in requested_reads:
                raise ValueError('Spec requests I7 read in barcode %s, but no index is found in the r1 name %s.' % (barcode, r1_name))

        # Get all barcode seqs
        for barcode in spec:
            bc_read = spec[barcode][BC_READ]
            bc_start = spec[barcode][BC_START] - 1
            bc_end = spec[barcode][BC_END]


            if bc_read == I5:
                read_seq = r1_name[i5_start:i5_end]
                if not read_lengths_checked:
                    _validate_barcode_read_pair(read_seq, bc_end)

                if reverse_i5:
                    bc_seq = reverse_complement(read_seq)[bc_start:bc_end]
                else:
                    bc_seq = read_seq[bc_start:bc_end]
 
            elif bc_read == I7:
                read_seq = r1_name[i7_start:i7_end]
                if not read_lengths_checked:
                    _validate_barcode_read_pair(read_seq, bc_end)
     
                bc_seq = read_seq[bc_start:bc_end]

            elif bc_read == R1:
                if not read_lengths_checked:
                    _validate_barcode_read_pair(r1_seq, bc_end)
                bc_seq = r1_seq[bc_start:bc_end]
                
            elif bc_read == R2:
                if not read_lengths_checked:
                    _validate_barcode_read_pair(r2_seq, bc_end)
                bc_seq = r2_seq[bc_start:bc_end]

            if BC_WHITELIST in spec[barcode]:
                bc_seq = spec[barcode][BC_CORRECTION_MAP].correct(bc_seq)

            barcodes_dict[barcode] = bc_seq

        read_lengths_checked = True
        yield barcodes_dict
from barcodeutils import *
import barcodeutils as bu
from nose.tools import assert_raises

def test_is_dna():
    dna = 'ATAGC'
    output = is_dna(dna)
    assert output == True

    dna = 'ANAGC'
    output= is_dna(dna)
    assert output == True

    dna = 'ATAGD'
    output = is_dna(dna)
    assert output == False

    dna = None
    assert_raises(ValueError, is_dna, dna)

def test_hamming_distance():
    a = 'ATACG'
    b = 'ATACT'
    output = hamming_distance(a, b)
    assert output == 1

    a = 'ATACG'
    b = 'AAACT'
    output = hamming_distance(a, b)
    assert output == 2

    a = 'ATACG'
    b = 'ATACG'
    output = hamming_distance(a, b)
    assert output == 0

    a = 'ATACG'
    b = 'ATACGA'
    assert_raises(ValueError, hamming_distance, a, b)
    assert_raises(ValueError, hamming_distance, b, a)
    assert output == 0

    a = 'ATACG'
    b = 'GAACT'
    output = hamming_distance(a, b, capdistance=2)
    assert output == 2

def test_reverse_complement():
    dna = 'ATCGC'
    expected = 'GCGAT'
    assert reverse_complement(dna) == expected

    dna = 'ATCGN'
    expected = 'NCGAT'
    assert reverse_complement(dna) == expected

def test_validate_barcode_spec_totally_off_inputs():

    spec = {}
    output, error = validate_barcode_spec(spec)
    assert output == False and error is not None

    spec = []
    output, error = validate_barcode_spec(spec)
    assert output == False and error is not None

    spec = ''
    output, error = validate_barcode_spec(spec)
    assert output == False and error is not None

def test_validate_barcode_spec_more_subtle():
    # Valid spec
    spec = {
        'cell_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: R1
        }
    }
    output, error = validate_barcode_spec(spec)
    print(error)
    assert output == True and error is None

    # Valid spec
    spec = {
        'cell_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: R1
        },
        'p5_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: I5
        }        
    }
    output, error = validate_barcode_spec(spec)
    assert output == True and error is None

    # Try to overwrite a reserved key in dict (other seqs the parser would output)
    spec = {
        'r1_seq': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: R1
        }       
    }
    output, error = validate_barcode_spec(spec)
    assert output == False and validate_barcode_spec(spec)

    # End is before start
    spec = {
        'r1_seq': {
            BC_START: 10,
            BC_END: 1,
            BC_READ: R1
        }       
    }
    output, error = validate_barcode_spec(spec)
    assert output == False and validate_barcode_spec(spec)

    # Invalid coord (1-based so anything 0 or below is not allowed)
    spec = {
        'r1_seq': {
            BC_START: 0,
            BC_END: 10,
            BC_READ: R1
        }       
    }
    output, error = validate_barcode_spec(spec)
    assert output == False and validate_barcode_spec(spec)

    # Missing required property
    spec = {
        'cell_barcode': {
            BC_START: 1,
            BC_END: 10
        },
        'p5_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: I5
        }        
    }
    output, error = validate_barcode_spec(spec)
    assert output == False and error is not None

    # Included sprurious property
    spec = {
        'cell_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: R1
        },
        'p5_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: I5,
            'INVALID': "ENTRY"
        }        
    }
    output, error = validate_barcode_spec(spec)
    assert output == False and error is not None

    # Included wrong property type
    spec = {
        'cell_barcode': {
            BC_START: '1',
            BC_END: 10,
            BC_READ: R1
        },
    }
    output, error = validate_barcode_spec(spec)
    assert output == False and error is not None

    # Another
    spec = {
        'cell_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: 1
        },
    }
    output, error = validate_barcode_spec(spec)
    assert output == False and error is not None

def test_validate_barcode_spec_whitelist():
    # Valid whitelist
    spec = {
        'cell_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: R1,
            BC_WHITELIST: ['ATACA', 'TATAC', 'ATACT']
        },
    }
    output, error = validate_barcode_spec(spec)
    assert output == True and error is None

    # Another valid whitelist
    spec = {
        'cell_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: R1,
            BC_WHITELIST: {'ATACA', 'TATAC', 'ATACT'}
        },
    }
    output, error = validate_barcode_spec(spec)
    assert output == True and error is None

    # Whitelist file
    spec = {
        'cell_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: R1,
            BC_WHITELIST: 'whitelist.txt'
        },
    }
    output, error = validate_barcode_spec(spec)
    assert output == True and error is None

    # File that doesn't exist
    spec = {
        'cell_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: R1,
            BC_WHITELIST: 'whereami.txt'
        },
    }
    output, error = validate_barcode_spec(spec)
    assert output == False and error is not None

def test_load_barcode_spec():
    # Load an example JSON file
    result = {
        'cell_barcode': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: R1,
            BC_WHITELIST: ['ATACA', 'TATAC', 'ATACT']
        },
    }
    output = load_barcode_spec('example_spec.json')
    assert output == result

    # Load a JSON file with a mistake (sets not allowed in Json)
    assert_raises(ValueError, load_barcode_spec, 'example_spec.bad.json')

def test_barcode_correcter():
    # Test correcter functionality
    correcter = BarcodeCorrecter(whitelist=['ATACA', 'TATAC', 'ATACT'], edit_distance=1)
    error_seq = 'TTACA'
    assert 'ATACA' == correcter.correct(error_seq)

    correcter = BarcodeCorrecter(whitelist=['ATACA', 'TATAC', 'ATACT'], edit_distance=1)
    error_seq = 'TTGCN'
    assert None == correcter.correct(error_seq)

    # Also test hamming distance
    correcter = BarcodeCorrecter(whitelist=['ATAA', 'AGAA', 'AACA'], edit_distance=1)
    assert [1] == correcter.get_min_hamming()
    
    correcter = BarcodeCorrecter(whitelist=['ATAA', 'AGAA', 'AACA'], edit_distance=1)
    assert [1,2] == correcter.get_min_hamming(n=2)

    ## Also request more pairs than exist, which should still work
    correcter = BarcodeCorrecter(whitelist=['ATAA', 'AGAA', 'AACA'], edit_distance=1)
    assert [1,2,2] == correcter.get_min_hamming(n=10)

    # Test barcode lengths
    correcter = BarcodeCorrecter(whitelist=['ATAA', 'AGAA', 'AACA'], edit_distance=1)
    assert 4 == correcter.get_barcode_length()


def test_get_index_coords():
    # Paired end
    r1_name = '...:N:0:TCGGATTCGG+CTCCATGGAG'
    i7_start = 8
    i7_end = 18
    i5_start = 19
    i5_end = 29
    output = bu._get_index_coords(r1_name)
    assert output == (i7_start, i7_end, i5_start, i5_end)

    ## Also test output seq to make sure is ok
    assert r1_name[i7_start:i7_end] == 'TCGGATTCGG'
    assert r1_name[i5_start:i5_end] == 'CTCCATGGAG'

    # Single ended
    r1_name = '...:N:0:TCGGATTCGG'
    i7_start = 8
    i7_end = 18
    i5_start = None
    i5_end = None
    output = bu._get_index_coords(r1_name)
    assert output == (i7_start, i7_end, i5_start, i5_end)

    ## Also test output seq to make sure is ok
    assert r1_name[i7_start:i7_end] == 'TCGGATTCGG'

    # Invalid (no indices)
    r1_name = '...:N:0'
    i7_start = None
    i7_end = None
    i5_start = None
    i5_end = None
    output = bu._get_index_coords(r1_name)
    assert output == (i7_start, i7_end, i5_start, i5_end)

def test_parse_fastq_barcodes():
    r1 = 'example.1.fastq.gz'
    r2 = 'example.2.fastq.gz'
    rt_barcodes = ["TTCTCGCATG", "TCCTACCAGT", "GCGTTGGAGC", "GATCTTACGC", "CTGATGGTCA", "CCGAGAATCC", "GCCGCAACGA", "TGAGTCTGGC", "TGCGGACCTA", "ACCTCGTTGA", "ACGGAGGCGG", "TAGATCTACT", "AATTAAGACT", "CCATTGCGTT", "TTATTCATTC", "ATCTCCGAAC"]
    i7_barcodes = ["TCGGATTCGG", "TCCGGCTTAT", "TCGCCGCCGG", "TTGGCAAGCC", "CAGCTAGCGG", "GTAGGATAAG", "AAGTAGCTCA", "TAAGTCCTGA", "GCGGCTGCGG", "ACCAGGCGCA", "CCGTATGATT", "TTGATTGGCG"]
    i5_barcodes = ["CTCCATCGAG", "TTGGTAGTCG", "GGCCGTCAAC", "CCTAGACGAG", "TCGTTAGAGC", "CGTTCTATCA", "CGGAATCTAA", "ATGACTGATC"]

    spec = {
        'umi': {
            BC_START: 1,
            BC_END: 8,
            BC_READ: R1
        },
        'cell_barcode': {
            BC_START: 9,
            BC_END: 18,
            BC_READ: R1,
            BC_WHITELIST: rt_barcodes
        },
        'p5': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: I5,
            BC_WHITELIST: i5_barcodes
        },
        'p7': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: I7,
            BC_WHITELIST: i7_barcodes
        }
    }

    reads = [x for x in parse_fastq_barcodes(r1, r2, spec=spec, reverse_i5=True, edit_distance=1)]
    assert reads[1]['p7'] == 'TCGGATTCGG'
    assert reads[1]['p5'] == 'CTCCATCGAG' # this read had a 1bp error to this barcode, should correct it
    assert reads[1]['umi'] == 'TATTTACC'
    assert reads[1]['cell_barcode'] == 'TGAGTCTGGC'

    # Accidentally don't reverse P5 when it had to be for this run
    reads = [x for x in parse_fastq_barcodes(r1, r2, spec=spec, reverse_i5=False, edit_distance=1)]
    assert reads[1]['p7'] == 'TCGGATTCGG'
    assert reads[1]['p5'] == None # this was a nextseq run, so don't expect to match the barcode
    assert reads[1]['umi'] == 'TATTTACC'
    assert reads[1]['cell_barcode'] == 'TGAGTCTGGC'

    # Exclude start and end coordinates
    spec = {
        'p5': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: I5,
            BC_WHITELIST: i5_barcodes
        },
        'p7': {
            BC_START: 1,
            BC_END: 10,
            BC_READ: I7,
            BC_WHITELIST: i7_barcodes
        }
    }

    reads = [x for x in parse_fastq_barcodes(r1, r2, spec=spec, reverse_i5=True, edit_distance=1)]
    assert reads[1]['p7'] == 'TCGGATTCGG'
    assert reads[1]['p5'] == 'CTCCATCGAG'

    # BC length doesn't match the whitelist length
    spec = {
        'p5': {
            BC_START: 1,
            BC_END: 8, # wrong
            BC_READ: I5,
            BC_WHITELIST: i5_barcodes
        },
        'p7': {
            BC_START: 1,
            BC_END: 8, # wrong
            BC_READ: I7,
            BC_WHITELIST: i7_barcodes
        }
    }

    def run_parse_fastq_barcodes(r1, r2, spec, reverse_i5, edit_distance):
        result = list(parse_fastq_barcodes(r1, r2, spec=spec, reverse_i5=reverse_i5, edit_distance=edit_distance))
        return result
    assert_raises(ValueError, run_parse_fastq_barcodes, r1, r2, spec, reverse_i5=True, edit_distance=1)

    # Requested region extend beyond read, but right length
    spec = {
        'p5': {
            BC_START: 3,
            BC_END: 12, # wrong
            BC_READ: I5,
            BC_WHITELIST: i5_barcodes
        },
        'p7': {
            BC_START: 3,
            BC_END: 12, # wrong
            BC_READ: I7,
            BC_WHITELIST: i7_barcodes
        }
    }
    assert_raises(ValueError, run_parse_fastq_barcodes, r1, r2, spec, reverse_i5=True, edit_distance=1)

def test_get_run_info():
    run_info = get_run_info('.')
    expected_output = {'flow_cell_id': 'HMHMGBGX7', 'date': '181011', 'instrument': 'NS500488', 'lanes': 4, 'r1_length': 18, 'p7_index_length': 10, 'p5_index_length': 10, 'r2_length': 52}
    assert run_info == expected_output

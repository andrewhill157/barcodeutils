# barcodeutils
In processing data from many new sequencing-based techs, you often have to parse out the sequences you care about, like indices, UMIs, barcodes in r1 or r2, etc. This also involves correcting some subset of these sequences to a known set of programmed sequences to allow for error correction.

`barcodeutils` automates each of these tasks with a simple interface.

## Example
For example, let's say you have an assay where you read i7, i5 reads as barcodes, the first 8bp of r1 as a unique molecular identifier (UMI), and the next 10 bp of r1 as another barcode.
```
import barcodeutils as bu

# Any known barcodes for error correction
my_r1_barcodes = ["TTCTCGCATG", "TCCTACCAGT", "GCGTTGGAGC", "GATCTTACGC", "CTGATGGTCA", "CCGAGAATCC", "GCCGCAACGA", "TGAGTCTGGC"]
my_i7_barcodes = ["TCGGATTCGG", "TCCGGCTTAT", "TCGCCGCCGG", "TTGGCAAGCC", "CAGCTAGCGG", "GTAGGATAAG", "AAGTAGCTCA", "TAAGTCCTGA", "GCGGCTGCGG", "ACCAGGCGCA", "CCGTATGATT", "TTGATTGGCG"]
my_i5_barcodes = ["CTCCATCGAG", "TTGGTAGTCG", "GGCCGTCAAC", "CCTAGACGAG", "TCGTTAGAGC", "CGTTCTATCA", "CGGAATCTAA", "ATGACTGATC"]

# Define specification
barcodes = {
    'i5': {
        'start': 1,
        'end': 10,
        'read': 'i5',
        'whitelist': my_i5_barcodes
    },
    'i7': {
        'start': 1,
        'end': 10,
        'read': 'i7',
        'whitelist': my_i7_barcodes
    },
    'umi': {
        'start': 1,
        'end': 8,
        'read': 'r1'
    },
    'r1_barcode': {
        'start': 9,
        'end': 18,
        'read': 'r1',
        'whitelist': my_r1_barcodes
    }
}

# Extract those corrected sequences from fastqs
fastq1 = 'myfastq.R1.fastq.gz'
fastq2 = 'myfastq.R2.fastq.gz'

for entry in bu.parse_fastq_barcodes(fastq1, fastq2, spec=barcodes, edit_distance=2):
    # you can now access any barcodes from spec
    umi_seq = entry['umi']
    r1_barcode_seq = entry['r1_barcode']

    # You can also access r1/r2 seq, name, qual automatically
    r1_seq = entry['r1_seq']
    r1_qual = entry['r1_qual']
    r1_name = entry['r1_name']
```

Any barcode defined with a whitelist will be corrected to the closest element within the whitelist within the specified edit distance. If there is no match or it would match to more than one barcode, then `None` is returned for that barcode seq.

## Storing Specifications and Barcode Whitelists
In the example above, the specification for the barcodes, etc. are hard-coded into the script, which may not be what you want. You can also store barcode whitelists in a text file with one 

You may put the whitelist sequences in a text file with one barcode per line:
```
TTCTCGCATG
TCCTACCAGT
GCGTTGGAGC
```

You may also put your barcode specification in a JSON file, for example:
```
{
    'i5': {
        'start': 1,
        'end': 10,
        'read': 'i5',
        'whitelist': 'my_i5_barcodes.txt'
    },
    'i7': {
        'start': 1,
        'end': 10,
        'read': 'i7',
        'whitelist': 'my_i7_barcodes.txt'
    },
    'umi': {
        'start': 1,
        'end': 8,
        'read': 'r1'
    },
    'r1_barcode': {
        'start': 9,
        'end': 18,
        'read': 'r1',
        'whitelist': 'my_r1_barcodes'
    }
}
```

Note that you may reference your whitelist text files, or you may include the sequences as a list directly in the spec, much like the original code example.

These JSON files can be loaded with:
```
barcode_spec = bu.load_barcode_spec('my_barcodes.json')
```

and any whitelist files will be loaded automatically.

## Dealing with Different Sequencers
A common problem is dealing with the fact that some sequencers, MiSeq/NextSeq sequencers for example, will produce i5 sequences that are a reverse complement of one another.

`parse_fastq_barcodes` has an additional argument `reverse_complement_i5` which will reverse complement the i5 read before extracting your specified i5 barcodes.

In case it is useful, `barcodeutils` also provides a helper function to help determine whether this needs to happen automatically taking a BCL directory (or a directory containing the `RunParameters.xml` file), or instrument type (NextSeq, MiSeq, NovaSeq, HiSeq4000, or HiSeq3000) as input. For any sequencers that behave like MiSeq, this will return False and for any sequencers that behave like NextSeq, this will return True. If you need for the orientation of your specification of whitelists, you may invert this as needed.

```
rc_i5 = reverse_complement_i5('NextSeq')
rc_i5 = reverse_complement_i5('/path/to/bcl')
```

## General Utilities
### Barcode Correction
If you just want to correct sequences to a known list of sequences, we also offer a standalone barcode correction tool.
```
my_r1_barcodes = ["TTCTCGCATG", "TCCTACCAGT", "GCGTTGGAGC", "GATCTTACGC", "CTGATGGTCA", "CCGAGAATCC", "GCCGCAACGA", "TGAGTCTGGC"]

r1_correcter = bu.BarcodeCorrecter(my_r1_barcodes, edit_distance=2)

corrected = r1_correcter.correct("CACCTTACGC") # GATCTTACGC

# Can also get other basic stats like the closest hamming distances or length
## TODO maybe make this spit back the pairs they correspond to
print(r1_correcter.get_min_hamming(n=3))

print(r1_correcter.get_barcode_length())
```

### Getting Run Info
If you have the BCL directory to the run and just want some basic info about it like read lengths, flow cell ID, date, instrument ID, lanes, you can get those returned as a dict like this:
```
bcl_directory = '/path/to/bcl/'
run_info = bu.get_run_info(bcl_directory)
```

This could also just be a directory containing the RunParameters.xml file from the BCL directory. The output will look something like the following, for example:
```
{
    'instrument_type': 'NextSeq',
    'lanes': 4,
    'instrument': 'NS500488',
    'p5_index_length': 0,
    'p7_index_length': 6,
    'flow_cell_id': 'HKFWJBGXX',
    'r2_length': 151,
    'r1_length': 151,
    'date': '151117'
}
```

### Reverse Complement DNA
```
seq = 'ATAGAGAC'
rc = bu.reverse_complement(seq)
```

## Installation
TODO I plan to put this on pip, but for now can install via the following or similar:
```
git clone git@github.com:andrewhill157/barcodeutils.git
cd barcodeutils
python setup.py install
```


"""
Make Datafile

This file will take in the twist dataset tables and splicing information and generate datafile.

The datafile will have the following columns in h5py format.
- ID
- Barcode
- Celltype
- full raw sequences, padded so everything is 900bp long
- Skipped Count
- Included Count

"""
import logging
import numpy as np
import sys
import time
import h5py
from constants import *
from utils import get_rev_comp

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)
# Length of the sequence.
DESIRED_LENGTH = 900
UPSTREAM_EXON = "ATCGCCACCATGGGTCGCTTTGACAGGGAGGTAGATATTGGAATTCCTGATGCTACAGGACGCTTAGAGATTCTTCAGATCCATACCAAGAACATGAAGCTGGCAGATGATGTGGACCTGGAACAG"
UPSTREAM_INTRON = "GTGAAGTGATGATGATGGCTGACCAGGCGTTACAGTGTCTCTAGGCAGTTGCTGGGAACTGGCTAGAGACATAAGGTTAAGATGTGAGGAGATGGGTTTTGATTTCTGGACAGGGGAAAGGAAGTAATCTGAGATTGAATCCAGGAAATGAAGCTTCGACACCGAGCTCG"
DOWNSTREAM_INTRON = "TGCAGCCATCTAAGTTTCATAAAGGCAGGGATTTTTGTCCATTTTAATTGATTTTTGTCAATTTCATTTAAAATTGTCATATAGTAAACACTCATTGTTTGTCCAAATGTCAAATGACTGGGTTTCCAGAACAGGTTGACACCTTTACTCCTGTTCTGTTACTGACCGGCCCTTCCCCTTCTGGGTCCTCCATTTCCTCATCTGGCACATGGACAGATGGCTGACCCATCCTGAAAGTCGCCTTGTTGTTCCTGCCCCAG"
DOWNSTREAM_EXON = "GTGCGGCAGCTGGTGCCTCGGAGGGACCAAGCTCTGACGGAGGAGCATGCCCGACAGCAGCACAATGAGAGGCTACGCAAGCAGTTTGGAGCCCAGGCCAATGTCATCGGGCCCTGGATCCAGACCAAG"

start_time = time.time()

assert sys.argv[1] in ["train", "test", "all"]
assert sys.argv[2] in ["0", "1", "all"]

if sys.argv[1] == "train":
    CHROM_GROUP = [
        "chr11",
        "chr13",
        "chr15",
        "chr17",
        "chr19",
        "chr21",
        "chr2",
        "chr4",
        "chr6",
        "chr8",
        "chr10",
        "chr12",
        "chr14",
        "chr16",
        "chr18",
        "chr20",
        "chr22",
        "chrX",
        "chrY",
    ]
elif sys.argv[1] == "test":
    CHROM_GROUP = ["chr1", "chr3", "chr5", "chr7", "chr9"]
else:
    CHROM_GROUP = [
        "chr1",
        "chr3",
        "chr5",
        "chr7",
        "chr9",
        "chr11",
        "chr13",
        "chr15",
        "chr17",
        "chr19",
        "chr21",
        "chr2",
        "chr4",
        "chr6",
        "chr8",
        "chr10",
        "chr12",
        "chr14",
        "chr16",
        "chr18",
        "chr20",
        "chr22",
        "chrX",
        "chrY",
    ]



CHROM = []
BARCODE = []
CELLTYPE = []
SEQ = []
SKIPPED_COUNT = []
INCLUDED_COUNT = []


seq_dict = {}
chrom_dict = {}

##############################################################################
# Read in the twist sequences.
##############################################################################
f = open(sequence_file)
lines = f.readlines()
header = lines.pop(0).strip().split(",")
for line in lines:
    (
        id,
        barcode,
        upstream_intron_seq,
        skipped_exon_seq,
        downstream_intron_seq,
        library_sequence,
        twist_sequence,
    ) = line.strip().split(",")
    # Split the id also. This ID is of format ENSG00000064932.16;SBNO2;chr19-1112188-1112301-1111995-1112067-1112401-1112537
    gene_id, gene_name, raw_id = id.split(";")
    # Get the chr number like "chr15" from chr19-4448317-4448415-4447549-4447625-4453457-4453522
    chr_number = raw_id.split("-")[0]

    if chr_number not in CHROM_GROUP:
        continue

    full_sequence = (
        UPSTREAM_EXON
        + UPSTREAM_INTRON
        + library_sequence
        + DOWNSTREAM_INTRON
        + DOWNSTREAM_EXON
    )
    # If the length is less than 900, pad with Ns on both sides.
    if len(full_sequence) < DESIRED_LENGTH:
        num_N_needed = DESIRED_LENGTH - len(full_sequence)
        right_N_needed = num_N_needed // 2
        left_N_needed = num_N_needed - right_N_needed
        full_sequence = "N" * left_N_needed + full_sequence + "N" * right_N_needed
    
    # Add to the dictionary.
    barcode = get_rev_comp(barcode)
    seq_dict[barcode] = full_sequence
    chrom_dict[barcode] = chr_number

f.close()

##############################################################################
# Read in the splicing results table.
##############################################################################
f = open(splicing_results_file)
lines = f.readlines()
header = lines.pop(0).strip().split(",")
n = 0
for line in lines:
    n += 1
    (
        _,
        barcode,
        unspliced_count,
        included_count,
        skipped_count,
        sample,
        celltype,
    ) = line.strip().split(",")
    # Strip the double quotes from barcode.
    barcode = barcode.strip('"')
    celltype = celltype.strip('"')
    if barcode in seq_dict:
        CHROM.append(chrom_dict[barcode])
        BARCODE.append(barcode)
        CELLTYPE.append(celltype)
        SEQ.append(seq_dict[barcode])
        skipped_count = skipped_count.strip('"')
        included_count = included_count.strip('"')
        skipped_count = int(skipped_count)
        included_count = int(included_count)
        SKIPPED_COUNT.append(skipped_count)
        INCLUDED_COUNT.append(included_count)
print(n)
f.close()

##############################################################################
# Save to h5py output.
##############################################################################
h5f = h5py.File(
    data_dir + "datafile" + "_" + sys.argv[1] + "_" + sys.argv[2] + ".h5", "w"
)

h5f.create_dataset("CHROM", data=np.asarray(CHROM).astype("|S"))
h5f.create_dataset("BARCODE", data=np.asarray(BARCODE).astype("|S"))
h5f.create_dataset("CELLTYPE", data=np.asarray(CELLTYPE).astype("|S"))
h5f.create_dataset("SEQ", data=np.asarray(SEQ).astype("|S"))
h5f.create_dataset("SKIPPED_COUNT", data=np.asarray(SKIPPED_COUNT))
h5f.create_dataset("INCLUDED_COUNT", data=np.asarray(INCLUDED_COUNT))

h5f.close()

print("--- %s seconds ---" % (time.time() - start_time))

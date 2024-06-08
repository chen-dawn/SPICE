"""
How to run the script.
python /broad/dawnccle/melange/process_fastq/missplicing/extract_fastq_for_K700E_shortlist.py 

"""
import gzip
import logging
import os
import pandas as pd
import sys
from Bio import SeqIO

def start_logger(outfile="RNA_MATCHING.log"):
    logging.basicConfig(
        datefmt="%H:%M:%S",
        handlers=[logging.FileHandler(outfile), logging.StreamHandler(sys.stdout)],
        format="%(asctime)s - %(message)s",
        level=logging.INFO,
    )

def extract_reads_from_fq(fq1_file, barcode_dict, output_dir):
    """
    Extract reads from a fastq file that match the barcodes in the provided dictionary.
    """
    num_reads = {barcode: 0 for barcode in barcode_dict}
    total_reads = 0
    reads_dict = {barcode: [] for barcode in barcode_dict}
    print(f"Extracting reads for {len(barcode_dict)} barcodes.")
    unzipped_file1 = gzip.open(fq1_file, "rt")
    for fq1 in SeqIO.parse(unzipped_file1, "fastq"):
        total_reads += 1
        if total_reads % 10000 == 0:
            print(f"Processed {total_reads} reads.")
        if total_reads >= 20000000:
            break

        header = fq1.id
        read_name = header.split()[0]
        id, cb, umi = read_name.split("_")
        if cb in barcode_dict:
            if num_reads[cb] < 100:
                num_reads[cb] += 1
                reads_dict[cb].append(fq1.format("fastq"))

    for barcode, reads in reads_dict.items():
        output_file = os.path.join(output_dir, f"{barcode_dict[barcode]}.fastq")
        with open(output_file, "w") as fo:
            fo.writelines(reads)
        logging.info(f"Extracted {num_reads[barcode]} reads for barcode: {barcode_dict[barcode]}")

shortlist_path = "/broad/dawnccle/processed_data/KMRC1_shortlist/high_KMRC1_shortlisted_elements.csv"
output_dir = "/broad/dawnccle/processed_data/KMRC1_shortlist/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

fastq_file = "/broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/merged_fastqs/KMRC1-rep1_R1_bc_extracted.fastq.gz"
# Read the shortlist.
shortlist = pd.read_csv(shortlist_path)
# Get the index and barcodeRevcomp columns.
shortlist = shortlist[["index", "barcodeRevcomp"]]
# Unique barcodes.
shortlist = shortlist.drop_duplicates()

# Prepare the barcode dictionary.
barcode_dict = {}
for i, row in shortlist.iterrows():
    index = row["index"]
    filename = "KMRC1-rep1_" + index.split(";")[1] + "_" + index.split(";")[2]
    barcode_dict[row["barcodeRevcomp"]] = filename

# Extract the reads for all barcodes in one pass.
extract_reads_from_fq(fastq_file, barcode_dict, output_dir)

# import os
# import mappy as mp
# from fnmatch import fnmatch
# from pathlib import Path
# import argparse
# import editdistance
# import sys
# from Bio import SeqIO
# import pandas as pd
# import gzip
# from collections import Counter
# import logging
# import edlib

# def start_logger(outfile="RNA_MATCHING.log"):
#     logging.basicConfig(
#         datefmt="%H:%M:%S",
#         handlers=[logging.FileHandler(outfile), logging.StreamHandler(sys.stdout)],
#         format="%(asctime)s - %(message)s",
#         level=logging.INFO,
#     )

# def extract_reads_from_fq(fq1_file, barcode, output_file):
#     """
#     Extract reads from a fastq file that match the barcode.
#     """
#     num_reads = 0
#     total_reads = 0
#     fo = open(output_file, "w")
#     unzipped_file1 = gzip.open(fq1_file, "rt")
#     for fq1 in SeqIO.parse(unzipped_file1, "fastq"):
#         total_reads += 1
#         if total_reads % 1000000 == 0:
#             logging.info(f"Processed {total_reads} reads.")
#         if total_reads > 10000000:
#             break
#         if num_reads >= 100:
#             break
#         header = fq1.id
#         read_name = header.split()[0]
#         id, cb, umi = read_name.split("_")
#         if cb == barcode:
#             num_reads += 1
#             print(num_reads)
#             fo.write(fq1.format("fastq"))
#     fo.close()
#     logging.info(f"Extracted {num_reads} reads for barcode: {barcode}")

# shortlist_path = "/broad/dawnccle/processed_data/K700E_shortlist/high_K700E_shortlisted_barcodes.csv"
# shortlist_path = "/broad/dawnccle/processed_data/K700E_shortlist/high_WT_shortlisted_elements.csv"
# shortlist_path = "/broad/dawnccle/processed_data/KMRC1_shortlist/high_KMRC1_shortlisted_elements.csv"
# output_dir = "/broad/dawnccle/processed_data/KMRC1_shortlist/"
# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)
# # fastq_file = "/broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/old_fastq/K562_K700E-H04_S262_R1_bc_extracted.fastq.gz"
# fastq_file = "/broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/merged_fastqs/KMRC1-rep1_R1_bc_extracted.fastq.gz"
# # Read the shortlist.
# shortlist = pd.read_csv(shortlist_path)
# # Get the index and barcodeRevcomp columns.
# shortlist = shortlist[["index", "barcodeRevcomp"]]
# # Unique barcodes.
# shortlist = shortlist.drop_duplicates()

# # Extract the reads for K700E.
# fastq_file = "/broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/old_fastq/K562_K700E-H04_S262_R1_bc_extracted.fastq.gz"
# fastq_file = "/broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/merged_fastqs/KMRC1-rep1_R1_bc_extracted.fastq.gz"
# for i, row in shortlist.iterrows():
#     print(f"Extracting reads for barcode: {row['index']}")
#     barcode = row["barcodeRevcomp"]
#     index = row["index"]
#     filename = "KMRC1-rep1_" + index.split(";")[1] + "_" + index.split(";")[2]
#     output_file = os.path.join(output_dir, f"{filename}.fastq")
#     extract_reads_from_fq(fastq_file, barcode, output_file)
#     logging.info(f"Extracted reads for barcode: {filename}")

# # Extract reads for WT. 
# # fastq_file = "/broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/old_fastq/K562_WT-H01_S259_R1_bc_extracted.fastq.gz"
# # for i, row in shortlist.iterrows():
# #     print(f"Extracting reads for barcode: {row['index']}")
# #     barcode = row["barcodeRevcomp"]
# #     index = row["index"]
# #     filename = "K562_WT-H01_" + index.split(";")[1] + "_" + index.split(";")[2]
# #     output_file = os.path.join(output_dir, f"{filename}.fastq")
# #     extract_reads_from_fq(fastq_file, barcode, output_file)
# #     logging.info(f"Extracted reads for barcode: {filename}")

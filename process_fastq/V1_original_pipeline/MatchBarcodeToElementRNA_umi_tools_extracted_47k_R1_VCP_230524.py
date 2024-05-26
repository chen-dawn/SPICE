"""
Matching sequenced RNA library to unspliced, included, and skipped.
Dawn Chen

Usage:

python MatchBarcodeToElementRNA_umi_tools_extracted_47k.py \
    -1 /broad/thechenlab/NextSeqOutput_Xerneas/230306_VL00297_150_AACG7YFM5/Data/Intensities/BaseCalls/A375-1-C04_S16_R1_bc_extracted.fastq.gz \
    -2 /broad/thechenlab/NextSeqOutput_Xerneas/230306_VL00297_150_AACG7YFM5/Data/Intensities/BaseCalls/A375-1-C04_S16_R2_bc_extracted.fastq.gz \
    -l /broad/thechenlab/Dawn/splicing/library_sequencing_47k/20230130_twist_library_v3_ID_barcode_ROUT.csv \
    -r /broad/thechenlab/Dawn/splicing/library_sequencing_47k/RNA_ref/WEAK_47k_reference.fasta \
    -o /broad/thechenlab/Dawn/splicing/library_sequencing_47k/230307_RNA_3celltype_try1


# The old code: 
python /broad/thechenlab/Dawn/splicing/library_sequencing_47k/MatchBarcodeToElementRNA_umi_tools_extracted_47k_R1_VCP_230524.py \
    -1 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/HEK-G06_S174_R1_bc_extracted.fastq.gz \
    -2 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/HEK-G06_S174_R2_bc_extracted.fastq.gz \
    -l /broad/thechenlab/Dawn/splicing/library_sequencing_47k/20230130_twist_library_v3_ID_barcode_ROUT.csv \
    -r /broad/thechenlab/Dawn/splicing/library_sequencing_47k/RNA_ref/WEAK_47k_reference.fasta \
    -o /broad/dawnccle/processed_data/celltype47_230524


python /broad/thechenlab/Dawn/splicing/library_sequencing_47k/MatchBarcodeToElementRNA_umi_tools_extracted_47k_R1_VCP_230524.py \
    -1 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/OC316-A12_S108_R1_bc_extracted.fastq.gz \
    -2 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/OC316-A12_S108_R2_bc_extracted.fastq.gz \
    -l /broad/thechenlab/Dawn/splicing/library_sequencing_47k/20230130_twist_library_v3_ID_barcode_ROUT.csv \
    -r /broad/thechenlab/Dawn/splicing/library_sequencing_47k/RNA_ref/WEAK_47k_reference.fasta \
    -o /broad/dawnccle/processed_data/celltype47_230524


To run locally: 
"""

import os
import mappy as mp
from fnmatch import fnmatch
from pathlib import Path
import argparse
import editdistance
import sys
from Bio import SeqIO

# import distance
import pandas as pd
import gzip
from collections import Counter
import logging


# logger = logging.getLogger('MatchBarcodeToElementRNA')
# logger.setLevel(logging.INFO)
# fh = logging.FileHandler('spam.log')
# fh.setLevel(logging.INFO)
# logger.addHandler(fh)

# Types of library.
WEAK = "WEAK"
MEDIUM = "MEDIUM"
STRONG = "STRONG"

def start_logger(outfile="RNA_MATCHING.log"):
    logging.basicConfig(
    datefmt="%H:%M:%S",
    handlers=[
        logging.FileHandler(outfile),
        logging.StreamHandler(sys.stdout)
    ],
    format="%(asctime)s - %(message)s",
    level=logging.INFO,
)

def parse_args():

    parser = argparse.ArgumentParser(
        description="Match RNA splicing sequencing to elements. "
    )

    parser.add_argument(
        "-1", "--fq1", type=str, required=True, help="Forward fastq.",
    )

    parser.add_argument(
        "-2", "--fq2", type=str, required=True, help="Reverse fastq. Assume that the umi and cell barcode are already in the header.",
    )

    parser.add_argument(
        "-n", "--sample_name", type=str, help="Name of the sample.", default="sample",
    )

    parser.add_argument(
        "-l",
        "--library_bc_path",
        type=str,
        required=True,
        help="The location of the barcode/library matching dictionary.",
    )

    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        help="Where the reference is located. Generate using the script MakeRNAElementDict.py",
        default="/broad/thechenlab/Dawn/splicing/library_sequencing_V2/RNA_ref/MEDIUM_reference.fasta",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help="Directory where output files will be saved.",
        default="/broad/thechenlab/Dawn/splicing/library_sequencing_V2/",
    )

    args = parser.parse_args()
    return args


def get_best_bc_from_dict(barcode, bc_lib_dict, max_dist=1):
    """
    Get closest barcode that's in the dictionary.
    """
    b_found = [
        b for b in bc_lib_dict.keys() if editdistance.eval(barcode, b) <= max_dist
    ]
    if len(b_found) == 0:
        return None
    if len(b_found) > 1:
        logging.info("Multiple matches detected: {}".format(len(b_found)))
    return b_found[0]


def get_identity_from_fq(fq1_file, fq2_file, guide_fasta_path, bc_lib_dict):
    aln = mp.Aligner(fn_idx_in=guide_fasta_path, preset="sr", best_n=2)
    total_reads = 0
    aligned_reads = 0
    bc_not_found = 0
    chimera_reads = 0
    chimera_reads_assigned = 0
    logging.info("Parsing fastq files.")

    # Create counters for UNSPLICED, INCLUDED, and SKIPPED.
    unspliced_counter_dict = {}
    included_counter_dict = {}
    skipped_counter_dict = {}
    umi_dedup_counter_dict = {}

    # Parse the fastq.
    unzipped_file1 = gzip.open(fq1_file,"rt")
    for fq1 in SeqIO.parse(unzipped_file1, "fastq"):
        is_aligned = False
        # print(fq1.seq)
        header = fq1.id
        read_name = header.split()[0]
        id, cb, umi = read_name.split("_")
        if cb not in umi_dedup_counter_dict:
            umi_dedup_counter_dict[cb] = [{}, {}, {}]


        total_reads += 1
        if total_reads % 10000 == 0:
            logging.info("Already processed reads: {}.".format(total_reads))

        # Extract read and aligned region.
        guide_bc = cb
        library_seq = str(fq1.seq)

        # We don't need to check anymore because the cell barcode is from the header whitelist. 
        if guide_bc not in bc_lib_dict:
            print(guide_bc, id, cb, umi)
            raise Exception("The barcode should exist in the whitelist.")

        element_id = bc_lib_dict[guide_bc]

        # Look to see if it's spliced or not?
        cb_umi = guide_bc + "_" + umi
        for al in aln.map(library_seq):
            if is_aligned:
                break
            if al.is_primary and al.mapq > 0:
                aligned_reads += 1
                is_aligned = True
                mapped_element = al.ctg
                # Now assign the element.
                if mapped_element == "EXON_SKIPPED":
                    if cb_umi not in skipped_counter_dict:
                        skipped_counter_dict[cb_umi] = 0
                    skipped_counter_dict[cb_umi] = skipped_counter_dict[cb_umi] + 1
                    umi_dedup_counter_dict[cb][2][umi] = 1

                elif mapped_element == "UNSPLICED":
                    if cb_umi not in unspliced_counter_dict:
                        unspliced_counter_dict[cb_umi] = 0
                    unspliced_counter_dict[cb_umi] = (
                        unspliced_counter_dict[cb_umi] + 1
                    )
                    umi_dedup_counter_dict[cb][0][umi] = 1
    
                elif mapped_element == element_id + "_INCLUDED":
                    if cb_umi not in included_counter_dict:
                        included_counter_dict[cb_umi] = 0
                    included_counter_dict[cb_umi] = (
                        included_counter_dict[cb_umi] + 1
                    )
                    umi_dedup_counter_dict[cb][1][umi] = 1
                else:
                    # Else we consider it as a chimera read.
                    # print("\n\nChimera read detected.")
                    # print(mapped_element)
                    # print(element_id)
                    # print(al)
                    chimera_reads += 1

    # Merge the elements into one dataframe.
    merged_df = pd.DataFrame(
        {
            "unspliced": pd.Series(unspliced_counter_dict),
            "included": pd.Series(included_counter_dict),
            "skipped": pd.Series(skipped_counter_dict),
        }
    )

    # Make the nested dictionary into a dataframe.
    for cb in umi_dedup_counter_dict:
        for i in range(3):
            umi_dedup_counter_dict[cb][i] = len(umi_dedup_counter_dict[cb][i])
    umi_dedup_df = pd.DataFrame(umi_dedup_counter_dict)

    # Add header names to the dataframe.
    umi_dedup_df = umi_dedup_df.transpose()
    # Add header names to the dataframe.
    umi_dedup_df.columns = ["unspliced", "included", "skipped"]

    logging.info("Total reads: {}".format(total_reads))
    logging.info("Total aligned reads: {}".format(aligned_reads))
    logging.info("Barcode not found reads: {}".format(bc_not_found))
    logging.info(
        "Percentage barcode not found reads: {}".format(
            bc_not_found / float(total_reads)
        )
    )
    logging.info("Chimera reads: {}".format(chimera_reads))
    logging.info(
        "Percentage chimera reads: {}".format(chimera_reads / float(total_reads))
    )
    logging.info("Chimera reads assigned: {}".format(chimera_reads_assigned))

    # Make a logs array
    stats_log = []
    stats_log.append(f"total_reads,{total_reads}\n")
    stats_log.append(f"total_aligned_reads,{aligned_reads}\n")
    stats_log.append(f"bc_not_found,{bc_not_found}\n")
    stats_log.append(f"perc_bc_not_found,{bc_not_found / float(total_reads)}\n")
    stats_log.append(f"chimera_reads,{chimera_reads}\n")
    stats_log.append(f"perc_chimera_reads,{chimera_reads / float(total_reads)}\n")
    return (merged_df, umi_dedup_df, stats_log)


def read_barcode_to_lib_table(table_path):
    d = {}
    f = open(table_path)
    lines = f.readlines()
    for line in lines:
        # Skip the header line.
        if "barcodeRevcomp" in line:
            continue
        id, barcode, barcode_rev_comp = line.strip().split(",")
        # print(id, barcode, barcode_rev_comp)
        # print(barcode_rev_comp)
        if barcode_rev_comp  == "GCTTCCGTGCTTGT":
            print("found it")
        d[barcode_rev_comp.strip('"')] = id.strip('"')
    return d


if __name__ == "__main__":
    args = parse_args()

    # sample_name = args.sample_name
    fq1_path = args.fq1
    fq2_path = args.fq2
    out_dir = args.output_dir

    # Create logging filename based on input file name.
    logging_filename = os.path.basename(fq1_path).replace(".fastq.gz", ".log")
    # Start logging.
    start_logger(logging_filename)


    # Create output dir if it doesn't exist.
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        logging.info("'{}' doesnt exist, create dir.".format(out_dir))

    logging.info("################# Starting on file: {}".format(Path(fq1_path)))
    
    # Read barcode to element matching.
    library_bc_path = args.library_bc_path
    logging.info("Reading barcode from library table: {}".format(library_bc_path))
    barcode_lib_dict = read_barcode_to_lib_table(library_bc_path)
    # print("barcode dict")
    # print(barcode_lib_dict)
    # Read guide reference dict.
    library_reference_dir = args.reference

    # This is the processed dataframe that's important.
    identity_df, umi_dedup_df, stats_log = get_identity_from_fq(
        fq1_path, fq2_path, library_reference_dir, barcode_lib_dict
    )

    # Write to output.
    the_basename = Path(fq1_path).stem
    the_basename = the_basename.replace("R1_bc_extracted.fastq", "")
    out_filename = the_basename + "splicing_counter.csv"
    outfile_path = os.path.join(out_dir, out_filename)
    logging.info("Writing output table to: {}".format(outfile_path))
    identity_df.to_csv(outfile_path)

    # Write umi dedup table.
    umi_dedup_outfile_path = os.path.join(out_dir, the_basename + "umi_dedup.csv")
    logging.info("Writing umi dedup table to: {}".format(umi_dedup_outfile_path))
    umi_dedup_df.to_csv(umi_dedup_outfile_path)

    # Write stats logs.
    stats_outfile_path = os.path.join(out_dir, the_basename + "stats_log.txt")
    f = open(stats_outfile_path, "w")
    for line in stats_log:
        f.write(line)
    f.close()


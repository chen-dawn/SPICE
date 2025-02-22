"""
This is now a new file that is going to try to find the microexons (<20 bp exons). We try our best! We are going to find the upstream and downstream constant exon sequences, and then just save the output. 

Dawn Chen

Usage:


HEK test file for debugging.
python /broad/dawnccle/melange/process_fastq/missplicing/Microexon_MatchBarcodeToElement_V1.py \
    -1 /broad/dawnccle/sequencing/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/merged_fastqs/HEK-rep1_R1_bc_extracted.fastq.gz \
    -l /broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv \
    -o /broad/dawnccle/processed_data/missplicing_debug
    
"""

import os
import mappy as mp
from fnmatch import fnmatch
from pathlib import Path
import argparse
import editdistance
import sys
from Bio import SeqIO
import pandas as pd
import gzip
from collections import Counter
import logging
import edlib

UPSTREAM_EXON_LAST_20 = "ATGATGTGGACCTGGAACAG"  # last 20 bp
DOWNSTREAM_EXON_FIRST_20 = "GTGCGGCAGCTGGTGCCTCG"  # first 20 bp
DOWNSTREAM_EXON_FIRST_20_3BP_OFFSET = (
    "CGGCAGCTGGTGCCTCGGAG"  # first 20 bp with 3 bp offset
)
DOWNSTREAM_EXON = "GTGCGGCAGCTGGTGCCTCGGAGGGACCAAGCTCTGACGGAGGAGCATGCCCGACAGCAGCACAATGAGAGGCTACGCAAGCAGTTTGGAGCCCAGGCCAATGTCATCGGGCCCTGGATCCAGACCAAG"
UNSPLICED_INTRON = "GTGAAGTGATGATGATGGCTGACCAGGCGTTACAGTGTCTCTAGGCAGTTGCTGGGAACTGGCTAGAGACATAAGGTTAAGATGTGAGGAGATGGGTTTTGATTTCTGGACAGGGGAAAGGAAGTAATCTGAGATTGAATCCAGGAAATGAAGCTTCGACACCGAGCTCG"
UPSTREAM_INTRON_LAST_40 = "GAGATTGAATCCAGGAAATGAAGCTTCGACACCGAGCTCG"  # 40 bp
DOWNSTREAM_INTRON = "TGCAGCCATCTAAGTTTCATAAAGGCAGGGATTTTTGTCCATTTTAATTGATTTTTGTCAATTTCATTTAAAATTGTCATATAGTAAACACTCATTGTTTGTCCAAATGTCAAATGACTGGGTTTCCAGAACAGGTTGACACCTTTACTCCTGTTCTGTTACTGACCGGCCCTTCCCCTTCTGGGTCCTCCATTTCCTCATCTGGCACATGGACAGATGGCTGACCCATCCTGAAAGTCGCCTTGTTGTTCCTGCCCCAG"


def start_logger(outfile="RNA_MATCHING.log"):
    logging.basicConfig(
        datefmt="%H:%M:%S",
        handlers=[logging.FileHandler(outfile), logging.StreamHandler(sys.stdout)],
        format="%(asctime)s - %(message)s",
        level=logging.INFO,
    )


def parse_args():

    parser = argparse.ArgumentParser(
        description="Match RNA splicing sequencing to elements. "
    )

    parser.add_argument(
        "-1",
        "--fq1",
        type=str,
        required=True,
        help="Forward fastq.",
    )

    parser.add_argument(
        "-2",
        "--fq2",
        type=str,
        help="Reverse fastq. Assume that the umi and cell barcode are already in the header.",
        default=None,
    )

    parser.add_argument(
        "-n",
        "--sample_name",
        type=str,
        help="Name of the sample.",
        default="sample",
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
        default="/broad/dawnccle/melange/data/guide_library_cleaned/WEAK_47k_reference_no_adapter_filtered.fasta",
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


def convert_csv_column2fa(csv_in, fa_out, id_col_name="id", fasta_col_name="fasta"):
    """
    Parses the specific column of a file into fasta format.
    Also returns the list of ids and sequences as a dictionary.
    """
    d = {}
    fo = open(fa_out, "w")
    f = open(csv_in)
    lines = f.readlines()
    id_number = 0
    column_number = 1
    header = lines.pop(0).strip().split(",")
    for i in range(len(header)):
        if header[i] == fasta_col_name:
            column_number = i
        if header[i] == id_col_name:
            id_number = i
    for line in lines:
        line_split = line.strip().split(",")
        fo.write(">{}\n{}\n".format(line_split[id_number], line_split[column_number]))
        d[line_split[id_number]] = line_split[column_number]
    f.close()
    fo.close()
    return d


def write_dict_to_fasta(guide_dict, fa_out):
    """
    Write a dictionary that has d[id] = seq into a fasta file.
    """
    fo = open(fa_out, "w")
    for id, seq in guide_dict.items():
        fo.write(">{}\n{}\n".format(id, seq))
    fo.close()


def read_guides(guide_input_csv, guide_fasta_out_path):
    guides_dict = convert_csv_column2fa(
        guide_input_csv,
        guide_fasta_out_path,
        id_col_name="id",
        fasta_col_name="librarySequence",
    )
    return guides_dict


def read_middle_exon(guide_input_csv, guide_middle_exon_out_path):
    middle_exon_dict = convert_csv_column2fa(
        guide_input_csv,
        guide_middle_exon_out_path,
        id_col_name="id",
        fasta_col_name="skippedExonSeq",
    )
    return middle_exon_dict


def read_upstream_intron(guide_input_csv, guide_upstream_intron_out_path):
    upstream_intron_dict = convert_csv_column2fa(
        guide_input_csv,
        guide_upstream_intron_out_path,
        id_col_name="id",
        fasta_col_name="upstreamIntronSeq",
    )
    return upstream_intron_dict


def read_downstream_intron(guide_input_csv, guide_downstream_intron_out_path):
    downstream_intron_dict = convert_csv_column2fa(
        guide_input_csv,
        guide_downstream_intron_out_path,
        id_col_name="id",
        fasta_col_name="downstreamIntronSeq",
    )
    return downstream_intron_dict


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


def find_substring_with_mismatches(dna_string, substring, max_mismatches):
    """
    Find substring with mismatches. Returns start and end index.
    """
    result = edlib.align(
        substring, dna_string, mode="HW", task="locations", k=max_mismatches
    )

    if result["locations"]:
        start, end = result["locations"][0]
        return start, end + 1  # edlib returns inclusive end index

    return None, None


def read_guide_fasta(guide_fasta_path):
    """
    Read the guide fasta file and return a dictionary.
    """
    guide_dict = {}
    with open(guide_fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                element_id = line.strip().replace(">", "")
                guide_dict[element_id] = ""
            else:
                guide_dict[element_id] = line.strip()
    return guide_dict


def add_cb_umi_to_dict(d, element, cb_umi):
    """
    This is a helper function to add the cell barcode and umi to the dictionary.
    """
    if element not in d:
        d[element] = {}
    if cb_umi not in d[element]:
        d[element][cb_umi] = 0
    d[element][cb_umi] += 1
    return d


def get_identity_from_fq(
    fq1_file,
    fq2_file,
    bc_lib_dict,
    guide_fasta_path,
    middle_exon_fasta_path,
    upstream_intron_fasta_path,
    downstream_intron_fasta_path,
):
    total_reads = 0
    aligned_reads = 0
    bc_not_found = 0

    element_cb_umi_dict = {}

    unzipped_file1 = gzip.open(fq1_file, "rt")
    for fq1 in SeqIO.parse(unzipped_file1, "fastq"):
        header = fq1.id
        read_name = header.split()[0]
        id, cb, umi = read_name.split("_")

        total_reads += 1
        
        if total_reads % 100000 == 0:
            logging.info(f"Processed {total_reads} reads.")

        guide_bc = cb
        cb_umi = cb + "_" + umi
        library_seq = str(fq1.seq)

        if guide_bc not in bc_lib_dict:
            bc_not_found += 1
            continue
        
        element_id = bc_lib_dict[guide_bc]
        
        # Find upstream exon boundary. Perfect match only. 
        start_upstream, end_upstream = find_substring_with_mismatches(
            library_seq, UPSTREAM_EXON_LAST_20, 0
        )
        if start_upstream is None:
            continue

        # Find downstream exon boundary. Perfect match only. 
        start_downstream, end_downstream = find_substring_with_mismatches(
            library_seq, DOWNSTREAM_EXON_FIRST_20, 0
        )
        if start_downstream is None:
            continue

        # Extract middle sequence between exon boundaries
        middle_seq = library_seq[end_upstream:start_downstream]
        if len(middle_seq) > 200:
            continue
        
        temp_ID = element_id + "_" + middle_seq
        
        aligned_reads += 1
        element_cb_umi_dict = add_cb_umi_to_dict(
            element_cb_umi_dict, temp_ID, cb_umi
        )
            
    # Convert to unnested dictionary with counts
    element_cb_umi_dict_unnested = {}
    for element in element_cb_umi_dict:
        element_cb_umi_dict_unnested[element] = len(element_cb_umi_dict[element])

    # Create dataframe and filter low counts
    umi_dedup_df = pd.DataFrame.from_dict(element_cb_umi_dict_unnested, orient="index")
    umi_dedup_df.columns = ["count"]
    umi_dedup_df.reset_index(level=0, inplace=True)
    umi_dedup_df = umi_dedup_df.sort_values(by="index")
    umi_dedup_df = umi_dedup_df[umi_dedup_df["count"] > 5]
    
    stats_log = []
    stats_log.append(f"total_reads,{total_reads}\n")
    stats_log.append(f"total_aligned_reads,{aligned_reads}\n")
    stats_log.append(f"bc_not_found,{bc_not_found}\n")

    return (umi_dedup_df, stats_log)


def read_barcode_to_lib_table(table_path):
    d = {}
    f = open(table_path)
    lines = f.readlines()
    for line in lines:
        if "barcodeRevcomp" in line:
            continue
        id, barcode, barcode_rev_comp = line.strip().split(",")
        d[barcode_rev_comp.strip('"')] = id.strip('"')
    return d

if __name__ == "__main__":
    guide_fasta_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_guides_filtered.fasta"
    middle_exon_fasta_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_skipped_exon_filtered.fasta"
    upstream_intron_fasta_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_upstream_intron_filtered.fasta"
    downstream_intron_fasta_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_downstream_intron_filtered.fasta"
    args = parse_args()

    fq1_path = args.fq1
    fq2_path = args.fq2
    out_dir = args.output_dir

    logging_filename = os.path.basename(fq1_path).replace(".fastq.gz", ".log")
    start_logger(logging_filename)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    library_bc_path = args.library_bc_path
    barcode_lib_dict = read_barcode_to_lib_table(library_bc_path)

    umi_dedup_df, stats_log = get_identity_from_fq(
        fq1_path,
        fq2_path,
        barcode_lib_dict,
        guide_fasta_path,
        middle_exon_fasta_path,
        upstream_intron_fasta_path,
        downstream_intron_fasta_path,
    )

    the_basename = Path(fq1_path).stem
    the_basename = the_basename.replace("R1_bc_extracted.fastq", "")

    umi_dedup_outfile_path = os.path.join(
        out_dir, the_basename + "microexon_umi_dedup_fine_grained_idx.csv"
    )
    umi_dedup_df.to_csv(umi_dedup_outfile_path, index=False)

    stats_outfile_path = os.path.join(
        out_dir, the_basename + "microexon_stats_log_fine_grained_idx.txt"
    )
    f = open(stats_outfile_path, "w")
    for line in stats_log:
        f.write(line)
    f.close()

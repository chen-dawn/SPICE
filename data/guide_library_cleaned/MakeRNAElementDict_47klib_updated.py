"""
Make RNA level element dictionaries. 
>WEAK
CAGGTGAAG	MAXENT: 6.66	
>MED
ATGGTATGT	MAXENT: 8.35	
>STRONG
ATGGTGAGT	MAXENT: 10.13	

Dawn Chen

"""

import os
import mappy as mp
from fnmatch import fnmatch
from pathlib import Path
# import distance
import pandas as pd
import gzip
from collections import Counter
import logging

logging.basicConfig(format="%(asctime)s - %(message)s", level=logging.INFO)

# Types of library.
WEAK = "weak"
MEDIUM = "medium"
STRONG = "strong"


WEAK_START = "GCTGGCAGATGATGTGGACCTGGAACAG"
WEAK_INTRON = "GTGAAGTGATGATGATGGCTGACCAGGCGTTACAGTGTCTCTAGGCAGTTGCTGGGAACTGGCTAGAGACATAAGGTTAAGATGTGAGGAGATGGGTTTTGATTTCTGGACAGGGGAAAGGAAGTAATCTGAGATTGAATCCAGGAAATG"
MEDIUM_START = "TGGGGCTGGACAGACTTTCTAATGAACCCAATG"  # There's 7Ns.
MEDIUM_INTRON = "GTATGTATCAGGCTAAAATAAAGCACTTTACTAATACCTAATTTGCTTCTAGGCAGTGGAAACTTCAAACCTCTGTGAGATAGAAAGCTGAGGTCAGATT"
STRONG_START = "AGGCAGACAAACCCATCAGCCATG"
STRONG_INTRON = "GTGAGTCTGCATCCTTTCCCCAGATGTGCCAATCATGGAGAGCCAGGCAGCAGCCACCACCATGCCCTGGAGTTGAGAGTAGAAGCTGTTGGAAAGATCATCTAACTGAGAAGAATTTTAATAGGGCATCAAAGATAAAGAATGCTGAGG"

DOWNSTREAM_EXON = "GTGCGGCAGCTGGTGCCTCGGAGGGACCAAGCTCTGACGGAGGAGCATGCCCGACAGCAGCACAATGAGAGGCTACGCAAGCAGTTTGGAGCCCAGGCCAATGTCATCGGGCCCTGGATCCAGACCAA"


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
    logging.info("Writing dict to fasta: {}".format(fa_out))
    fo = open(fa_out, "w")
    for id, seq in guide_dict.items():
        fo.write(">{}\n{}\n".format(id, seq))
    fo.close()


def read_guides(guide_input_csv, guide_fasta_out_path):
    guides_dict = convert_csv_column2fa(
        guide_input_csv,
        guide_fasta_out_path,
        id_col_name="ID",
        fasta_col_name="librarySequence",
    )
    return guides_dict


def read_skipped_exon(guide_input_csv, guide_skipped_exon_out_path):
    skipped_exon_dict = convert_csv_column2fa(
        guide_input_csv,
        guide_skipped_exon_out_path,
        id_col_name="ID",
        fasta_col_name="skippedExonSeq",
    )
    return skipped_exon_dict

def read_upstream_intron(guide_input_csv, guide_upstream_intron_out_path):
    upstream_intron_dict = convert_csv_column2fa(
        guide_input_csv,
        guide_upstream_intron_out_path,
        id_col_name="ID",
        fasta_col_name="upstreamIntronSeq",
    )
    return upstream_intron_dict


def read_downstream_intron(guide_input_csv, guide_downstream_intron_out_path):
    downstream_intron_dict = convert_csv_column2fa(
        guide_input_csv,
        guide_downstream_intron_out_path,
        id_col_name="ID",
        fasta_col_name="downstreamIntronSeq",
    )
    return downstream_intron_dict


def write_guide_bc_dict(aln_dict, outpath):
    logging.info("Writing barcode dict to file: {}".format(outpath))
    raw_f = open(outpath, "w")
    raw_f.write("id,barcode,count\n")
    for id in aln_dict:
        bcs = aln_dict[id]
        for bc in bcs:
            id_write = id.decode("UTF-8")
            bc_write = bc[0].decode("UTF-8")
            raw_f.write(f"{id_write},{bc_write},{bc[1]}\n")
    raw_f.close()


def make_RNA_splicing_reference(guide_dict, skipped_exon_dict, dict_type=MEDIUM):
    """"
    Make a RNA based reference. Going to make: 
    * Exon included version.
    * Unspliced version.
    """
    if dict_type == WEAK:
        upstream_exon = WEAK_START
        upstream_intron = WEAK_INTRON
    elif dict_type == MEDIUM:
        upstream_exon = MEDIUM_START
        upstream_intron = MEDIUM_INTRON
    elif dict_type == STRONG:
        upstream_exon = STRONG_START
        upstream_intron = STRONG_INTRON

    ref = {}
    skipped_version = upstream_exon + DOWNSTREAM_EXON
    ref["EXON_SKIPPED"] = skipped_version
    full_unspliced_seq = upstream_exon + upstream_intron
    ref["UNSPLICED"] = full_unspliced_seq

    for id in guide_dict:
        library_seq = guide_dict[id]
        skipped_exon_seq = skipped_exon_dict[id]

        full_included_seq = skipped_exon_seq + DOWNSTREAM_EXON

        ref[id + "_INCLUDED"] = full_included_seq[:180]

    return ref

if __name__ == "__main__":
    # Read guide dict and make reference
    reference_dir = "/broad/dawnccle/melange/data/guide_library_cleaned/"
    guide_input_csv = "/broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_filtered.csv"
    guide_fasta_out_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_guides_filtered.fasta"
    skipped_exon_fasta_out_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_skipped_exon_filtered.fasta"
    upstream_intron_fasta_out_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_upstream_intron_filtered.fasta"
    downstream_intron_fasta_out_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_downstream_intron_filtered.fasta"
    
    guide_dict = read_guides(guide_input_csv, guide_fasta_out_path)
    skipped_exon_dict = read_skipped_exon(guide_input_csv, skipped_exon_fasta_out_path)
    upstream_intron_dict = read_upstream_intron(guide_input_csv, upstream_intron_fasta_out_path)
    downstream_intron_dict = read_downstream_intron(guide_input_csv, downstream_intron_fasta_out_path)
    
    # Create reference directory and outpath directory.
    if not os.path.exists(reference_dir):
        os.makedirs(reference_dir)
        
    # Weak.
    library_reference = make_RNA_splicing_reference(
        guide_dict, skipped_exon_dict, WEAK
    )
    library_reference_dir = os.path.join(reference_dir, "WEAK_47k_reference_no_adapter_filtered.fasta")
    write_dict_to_fasta(library_reference, library_reference_dir)
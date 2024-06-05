"""
Matching sequenced RNA library to fine grained unspliced, included, and skipped.
Dawn Chen

Usage:

python /broad/dawnccle/melange/process_fastq/missplicing/MatchBarcodeToElementRNA_umi_tools_extracted_Novaseq230524_missplicing.py \
    -1 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/old_fastq/K562_K700E-H04_S262_R1_bc_extracted.fastq.gz \
    -2 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/old_fastq/K562_K700E-H04_S262_R2_bc_extracted.fastq.gz \
    -l /broad/dawnccle/melange/data/guide_library/20230130_twist_library_v3_ID_barcode_ROUT.csv \
    -r /broad/dawnccle/melange/data/guide_library/WEAK_47k_reference_no_adapter.fasta \
    -o /broad/dawnccle/processed_data/missplicing_test

python /broad/dawnccle/melange/process_fastq/missplicing/MatchBarcodeToElementRNA_umi_tools_extracted_Novaseq230524_missplicing.py \
    -1 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/test/DLST_K700E_R1_clean.fastq.gz \
    -2 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/old_fastq/K562_K700E-H04_S262_R2_bc_extracted.fastq.gz \
    -l /broad/dawnccle/melange/data/guide_library/20230130_twist_library_v3_ID_barcode_ROUT.csv \
    -r /broad/dawnccle/melange/data/guide_library_cleaned/WEAK_47k_reference_no_adapter_filtered.fasta \
    -o /broad/dawnccle/processed_data/missplicing_K700E


python /broad/dawnccle/melange/process_fastq/missplicing/MatchBarcodeToElementRNA_umi_tools_extracted_Novaseq230524_missplicing.py \
    -1 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/merged_fastqs/HEK-rep1_R1_bc_extracted.fastq.gz \
    -2 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/old_fastq/K562_K700E-H04_S262_R2_bc_extracted.fastq.gz \
    -l /broad/dawnccle/melange/data/guide_library/20230130_twist_library_v3_ID_barcode_ROUT.csv \
    -r /broad/dawnccle/melange/data/guide_library/WEAK_47k_reference_no_adapter.fasta \
    -o /broad/dawnccle/processed_data/missplicing_test
    

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

# Types of library.
WEAK = "WEAK"
MEDIUM = "MEDIUM"
STRONG = "STRONG"

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
        default=None
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
    # Perform the alignment
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
    spliced_reference,
    bc_lib_dict,
    guide_fasta_path,
    middle_exon_fasta_path,
    upstream_intron_fasta_path,
    downstream_intron_fasta_path,
):
    middle_exon_fasta = read_guide_fasta(middle_exon_fasta_path)
    upstream_intron_fasta = read_guide_fasta(upstream_intron_fasta_path)
    downstream_intron_fasta = read_guide_fasta(downstream_intron_fasta_path)
    guide_fasta = read_guide_fasta(guide_fasta_path)
    
    # Create an aligner.
    """
    Class mappy.Alignment
    This class describes an alignment. An object of this class has the following properties:

    ctg: name of the reference sequence the query is mapped to
    ctg_len: total length of the reference sequence
    r_st and r_en: start and end positions on the reference
    q_st and q_en: start and end positions on the query
    strand: +1 if on the forward strand; -1 if on the reverse strand
    mapq: mapping quality
    blen: length of the alignment, including both alignment matches and gaps but excluding ambiguous bases.
    mlen: length of the matching bases in the alignment, excluding ambiguous base matches.
    NM: number of mismatches, gaps and ambiguous positions in the alignment
    trans_strand: transcript strand. +1 if on the forward strand; -1 if on the reverse strand; 0 if unknown
    is_primary: if the alignment is primary (typically the best and the first to generate)
    read_num: read number that the alignment corresponds to; 1 for the first read and 2 for the second read
    cigar_str: CIGAR string
    cigar: CIGAR returned as an array of shape (n_cigar,2). The two numbers give the length and the operator of each CIGAR operation.
    MD: the MD tag as in the SAM format. It is an empty string unless the MD argument is applied when calling mappy.Aligner.map().
    cs: the cs tag.
    """
    aln = mp.Aligner(
        fn_idx_in=guide_fasta_path, preset="sr", best_n=3,k=9,w=1,min_chain_score=20)
    total_reads = 0
    aligned_reads = 0
    bc_not_found = 0
    included_reads_pre_chimera = 0
    chimera_reads = 0
    chimera_reads_assigned = 0

    perfect_skipped = 0
    imperfect_skipped = 0
    unspliced_counter = 0
    minor_splice_counter = 0
    no_primary_alignment = 0
    no_primary_alignment_but_mapped = 0
    multiple_primary_alignment = 0
    included_mapped_reads_but_low_quality = 0
    no_downstream_constant_exon = 0
    has_downstream_exon_offset = 0
    no_mapped_back_exon = 0
    perfect_middle_exon = 0
    everything_passed_with_middle_exon = 0
    logging.info("Parsing fastq files.")

    # Create a new dictionary to store the counts with splice information.
    element_cb_umi_dict = {}
    
    # Hits to merge. Seems to be some that has no primary alignment but are still mapped perfectly. 
    hits_to_merge_dict = {}

    # Parse the fastq. We only need to look at R1 since we already have the barcode and UMI.
    unzipped_file1 = gzip.open(fq1_file, "rt")
    # unzipped_file2 = gzip.open(fq2_file,"rt")
    for fq1 in SeqIO.parse(unzipped_file1, "fastq"):
        # print(fq1.seq)
        header = fq1.id
        read_name = header.split()[0]
        id, cb, umi = read_name.split("_")

        total_reads += 1
        if total_reads % 10000 == 0:
            logging.info("Already processed reads: {}.".format(total_reads))

        if total_reads >= 100000:
            break
        # Extract read and aligned region.
        guide_bc = cb
        cb_umi = cb + "_" + umi
        library_seq = str(fq1.seq)

        # We don't need to check anymore because the cell barcode is from the header whitelist.
        if guide_bc not in bc_lib_dict:
            bc_not_found += 1
            # print(guide_bc, id, cb, umi)
            # raise Exception("The barcode should exist in the whitelist.")
            continue
        
        element_id = bc_lib_dict[guide_bc]

        # Find the start and end of the "upstream exon" and "downstream exon" in the library sequence.
        # We will use this to determine if the read is spliced or not.
        start_upstream, end_upstream = find_substring_with_mismatches(
            library_seq, UPSTREAM_EXON_LAST_20, 2
        )
        start_downstream, end_downstream = find_substring_with_mismatches(
            library_seq, DOWNSTREAM_EXON_FIRST_20, 2
        )

        # If the upstream exon is not found, then we can't do anything.
        # We are going to skip this, there is 4 out of 10,000 reads that don't have this.
        if start_upstream is None or end_upstream is None:
            continue

        # Trim the library sequence to exclude the upstream exon now that we have it mapped.
        library_seq_upstream_exon_trimmed = library_seq[end_upstream:]
        # Get the first 20 or whatever is left.
        if len(library_seq_upstream_exon_trimmed) < 20:
            library_seq_first_20 = library_seq_upstream_exon_trimmed
        else:
            library_seq_first_20 = library_seq_upstream_exon_trimmed[:20]
            
        # Also get library_seq_first_40
        if len(library_seq_upstream_exon_trimmed) < 40:
            library_seq_first_40 = library_seq_upstream_exon_trimmed
        else:
            library_seq_first_40 = library_seq_upstream_exon_trimmed[:40]

        ########################################
        # Check out unspliced reads.
        ########################################
        # Try to map to the unspliced reference.
        unsplice_start, unsplice_end = find_substring_with_mismatches(
            UNSPLICED_INTRON, library_seq_first_20, 2
        )
        if unsplice_start is not None:
            aligned_reads += 1
            temp_ID = element_id + "_UNSPLICED_" + str(unsplice_start)
            element_cb_umi_dict = add_cb_umi_to_dict(
                element_cb_umi_dict, temp_ID, cb_umi
            )
            if unsplice_start != 0:
                minor_splice_counter += 1
            else:
                unspliced_counter += 1
            continue

        ########################################
        # Check out exon skipped reads.
        ########################################
        start_downstream_ref, end_downstream_ref = find_substring_with_mismatches(
            UPSTREAM_EXON_LAST_20 + DOWNSTREAM_EXON, library_seq_first_20, 2
        )
        if start_downstream_ref is not None and end_downstream_ref is not None:
            if abs(start_downstream_ref - len(UPSTREAM_EXON_LAST_20)) < 5:
                aligned_reads += 1
                difference = start_downstream_ref - len(UPSTREAM_EXON_LAST_20)
                temp_ID = element_id + "_SKIPPED_" + str(difference)
                # print(temp_ID)
                element_cb_umi_dict = add_cb_umi_to_dict(
                    element_cb_umi_dict, temp_ID, cb_umi
                )
                if difference == 0:
                    perfect_skipped += 1
                else:
                    imperfect_skipped += 1
                continue

        ########################################
        # Check out included reads.
        ########################################
        # Now we map.
        # First check how many primary alignments there are.
        num_primary = 0
        # We will just map the first 40 bp because if it's too long it stops working?
        # for al in aln.map(library_seq_upstream_exon_trimmed):
        
        # print(library_seq_first_40)
        for al in aln.map(library_seq_first_40):
            if al.is_primary and al.mapq > 0:
                print(al)
                num_primary += 1
                mapped_element = al.ctg
                
        # If there are no primary alignments, we check the alignment and give it to the alphabetically first element.
        if num_primary == 0:
            no_primary_alignment += 1
            continue
        
        if num_primary > 1:
            multiple_primary_alignment += 1
            # Print the multiple primary alignments.
            print("##############################################")
            for al in aln.map(library_seq_first_40):
                if al.is_primary and al.mapq > 0:
                    print(al)
            continue
        

        # Now we can process the primary alignment.
        aligned_reads += 1
        # Check if the mapping is consistent with the barcode. We will just keep track and not do anything right now.
        included_reads_pre_chimera += 1
        if mapped_element != element_id:
            chimera_reads += 1
        mapped_guide_seq = guide_fasta[mapped_element]
        mapped_upstream_intron_seq = upstream_intron_fasta[mapped_element]
        mapped_downstream_intron_seq = downstream_intron_fasta[mapped_element]
        mapped_middle_exon_seq = middle_exon_fasta[mapped_element]

        # Get the start and end of library_seq_first_20 in the mapped guide sequence.
        start_guide, end_guide = find_substring_with_mismatches(
            UPSTREAM_INTRON_LAST_40 + mapped_guide_seq, library_seq_first_20, 3
        )

        # TODO(dawnxchen): We are going to skip this for now, but actually these reads are super interesting. They seem to be spliced but not in the way that you expect.
        # Some of them are "correct" but it's just that there are more mismatches than 2. Not sure if it's DNA synthesis error.
        if start_guide is None:
            included_mapped_reads_but_low_quality += 1
            continue

        start_guide_pos_adjusted = (
            start_guide - len(mapped_upstream_intron_seq) - len(UPSTREAM_INTRON_LAST_40)
        )

        # Now we try to find the end of the downstream exon.
        last_20_bp_before_constant = None
        downstream_exon_offset = 0
        start_downstream, end_downstream = find_substring_with_mismatches(
            library_seq, DOWNSTREAM_EXON_FIRST_20, 2
        )
        if start_downstream is None:
            # Try to find the seq with 3 bp offset.
            start_downstream_w_3bp_offset, end_downstream_w_3bp_offset = (
                find_substring_with_mismatches(
                    library_seq, DOWNSTREAM_EXON_FIRST_20_3BP_OFFSET, 2
                )
            )

            if start_downstream_w_3bp_offset is not None:
                has_downstream_exon_offset += 1
                last_20_bp_before_constant = library_seq[
                    start_downstream_w_3bp_offset - 20 : start_downstream_w_3bp_offset
                ]
                downstream_exon_offset = 3
        else:
            last_20_bp_before_constant = library_seq[
                start_downstream - 20 : start_downstream
            ]

        # If we can't find the downstream exon, then we can't do anything.
        if last_20_bp_before_constant is None:
            no_downstream_constant_exon += 1
            continue

        # Now we try to map the last 20 bp before the constant exon to the downstream intron.
        # Increasing the error tolerance here since the read is worse quality around here now.
        start_downstream_exon, end_downstream_exon = find_substring_with_mismatches(
            UPSTREAM_INTRON_LAST_40 + mapped_guide_seq + DOWNSTREAM_INTRON,
            last_20_bp_before_constant,
            3,
        )

        if start_downstream_exon is None:
            no_mapped_back_exon += 1
            continue

        # Get the coordinate of the end of the downstream exon.
        end_guide_pos_adjusted = (
            end_downstream_exon
            - len(UPSTREAM_INTRON_LAST_40)
            - len(mapped_upstream_intron_seq)
            - len(mapped_middle_exon_seq)
        )
        # print(start_guide_pos_adjusted, end_guide_pos_adjusted)

        # Construct the name of the element.
        temp_ID = f"{mapped_element}_INCLUDED_{start_guide_pos_adjusted}:{end_guide_pos_adjusted}:{downstream_exon_offset}"
        element_cb_umi_dict = add_cb_umi_to_dict(element_cb_umi_dict, temp_ID, cb_umi)
        if start_guide_pos_adjusted == 0 and end_guide_pos_adjusted == 0:
            perfect_middle_exon += 1
        everything_passed_with_middle_exon += 1

    # print(element_cb_umi_dict)
    # Make a new dictionary that's the unnested element_cb_umi_dict.
    element_cb_umi_dict_unnested = {}
    for element in element_cb_umi_dict:
        element_cb_umi_dict_unnested[element] = len(element_cb_umi_dict[element])

    # Convert element_cb_umi_dict_unnested to a dataframe.
    umi_dedup_df = pd.DataFrame.from_dict(element_cb_umi_dict_unnested, orient="index")
    # Set the column name.
    umi_dedup_df.columns = ["count"]
    # Set the index as a column.
    umi_dedup_df.reset_index(level=0, inplace=True)
    # Sort by index.
    umi_dedup_df = umi_dedup_df.sort_values(by="index")
    # print(merged_df)

    # Convert element_cb_umi_dict to a dataframe. This is a nested dictionary of dictionaries so we unnest it.
    names_list = []
    cb_umis_list = []
    values_list = []
    for element in element_cb_umi_dict:
        for cb_umi in element_cb_umi_dict[element]:
            names_list.append(element)
            cb_umis_list.append(cb_umi)
            values_list.append(element_cb_umi_dict[element][cb_umi])
    # Create a data frame.
    element_cb_umi_df = pd.DataFrame(
        {
            "element": names_list,
            "cb_umi": cb_umis_list,
            "count": values_list,
        }
    )
    # Sort by element and cb_umi.
    element_cb_umi_df = element_cb_umi_df.sort_values(by=["element", "cb_umi"])
    
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
        "Percentage chimera reads: {}".format(
            chimera_reads / float(included_reads_pre_chimera)
        )
    )
    logging.info("Perfect skipped: {}".format(perfect_skipped))
    logging.info("Imperfect skipped: {}".format(imperfect_skipped))
    logging.info("Unspliced: {}".format(unspliced_counter))
    logging.info("Minor splice: {}".format(minor_splice_counter))
    logging.info("No primary alignment: {}".format(no_primary_alignment))
    logging.info("No primary alignment but mapped: {}".format(no_primary_alignment_but_mapped))
    logging.info("Multiple primary alignment: {}".format(multiple_primary_alignment))
    logging.info(
        "Included mapped reads but low quality: {}".format(
            included_mapped_reads_but_low_quality
        )
    )
    logging.info("No downstream constant exon: {}".format(no_downstream_constant_exon))
    logging.info("Has downstream exon offset: {}".format(has_downstream_exon_offset))
    logging.info("No mapped back exon: {}".format(no_mapped_back_exon))
    logging.info("Perfect middle exon: {}".format(perfect_middle_exon))
    logging.info(
        "Everything passed with middle exon: {}".format(
            everything_passed_with_middle_exon
        )
    )

    # Make a logs array
    stats_log = []
    stats_log.append(f"total_reads,{total_reads}\n")
    stats_log.append(f"total_aligned_reads,{aligned_reads}\n")
    stats_log.append(f"bc_not_found,{bc_not_found}\n")
    stats_log.append(f"perc_bc_not_found,{bc_not_found / float(total_reads)}\n")
    stats_log.append(f"chimera_reads,{chimera_reads}\n")
    stats_log.append(
        f"perc_chimera_reads,{chimera_reads / float(included_reads_pre_chimera)}\n"
    )
    stats_log.append(f"perfect_skipped,{perfect_skipped}\n")
    stats_log.append(f"imperfect_skipped,{imperfect_skipped}\n")
    stats_log.append(f"unspliced_counter,{unspliced_counter}\n")
    stats_log.append(f"minor_splice_counter,{minor_splice_counter}\n")
    stats_log.append(f"no_primary_alignment,{no_primary_alignment}\n")
    stats_log.append(f"no_primary_alignment_but_mapped,{no_primary_alignment_but_mapped}\n")
    stats_log.append(f"multiple_primary_alignment,{multiple_primary_alignment}\n")
    stats_log.append(
        f"included_mapped_reads_but_low_quality,{included_mapped_reads_but_low_quality}\n"
    )
    stats_log.append(f"no_downstream_constant_exon,{no_downstream_constant_exon}\n")
    stats_log.append(f"has_downstream_exon_offset,{has_downstream_exon_offset}\n")
    stats_log.append(f"no_mapped_back_exon,{no_mapped_back_exon}\n")
    stats_log.append(f"perfect_middle_exon,{perfect_middle_exon}\n")
    stats_log.append(
        f"everything_passed_with_middle_exon,{everything_passed_with_middle_exon}\n"
    )

    return (element_cb_umi_df, umi_dedup_df, stats_log)


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
        # if barcode_rev_comp == "GCTTCCGTGCTTGT":
        #     print("found it")
        d[barcode_rev_comp.strip('"')] = id.strip('"')
    return d


if __name__ == "__main__":
    guide_fasta_path = (
        "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_guides_filtered.fasta"
    )
    middle_exon_fasta_path = (
        "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_skipped_exon_filtered.fasta"
    )
    upstream_intron_fasta_path = (
        "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_upstream_intron_filtered.fasta"
    )
    downstream_intron_fasta_path = (
        "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_downstream_intron_filtered.fasta"
    )
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
        fq1_path,
        fq2_path,
        library_reference_dir,
        barcode_lib_dict,
        guide_fasta_path,
        middle_exon_fasta_path,
        upstream_intron_fasta_path,
        downstream_intron_fasta_path,
    )

    # Write to output.
    the_basename = Path(fq1_path).stem
    the_basename = the_basename.replace("R1_bc_extracted.fastq", "")
    out_filename = the_basename + "element_cb_umi_count_fine_grained_idx.csv"
    outfile_path = os.path.join(out_dir, out_filename)
    logging.info("Writing output table to: {}".format(outfile_path))
    identity_df.to_csv(outfile_path, index=False)

    # Write umi dedup table.
    umi_dedup_outfile_path = os.path.join(
        out_dir, the_basename + "umi_dedup_fine_grained_idx.csv"
    )
    logging.info("Writing umi dedup table to: {}".format(umi_dedup_outfile_path))
    umi_dedup_df.to_csv(umi_dedup_outfile_path, index=False)

    # Write stats logs.
    stats_outfile_path = os.path.join(
        out_dir, the_basename + "stats_log_fine_grained_idx.txt"
    )
    f = open(stats_outfile_path, "w")
    for line in stats_log:
        f.write(line)
    f.close()

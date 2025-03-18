"""
Script for analyzing RNA sequencing data from a splicing library. 
Takes fastq files containing RNA reads and matches them to a reference library of splicing elements. 

The code:
1. Extracts cell barcodes and UMIs from read headers
2. Maps reads to reference sequences to identify splicing patterns:
    - Unspliced reads that retain the intron
    - Skipped exon reads that splice out the middle exon
    - Included reads that retain the middle exon. 
    This works by aligning the "included region" of the read and then mapping 
3. Tracks statistics like chimeric reads and splicing efficiency
4. Outputs count tables of splicing patterns per cell/UMI

The unspliced read has one number that specifies the position of the intron where splicing occurs.
The skipped exon read has one numbers:
    - the coordinate on the constant downstream exon where splicing occurs.
The included read has three numbers:
    - the coordinate on the constant downstream exon where splicing occurs.
    - the start coordinate of the middle exon relative to the reference middle exon coordinate.
    - the end coordinate of the middle exon relative to the reference middle exon coordinate.

Then the 4 element is the insert length. 

Dawn Chen

Usage:

HEK test file for debugging.
python /broad/dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/MatchBarcodeToElementRNA_20250220.py \
    -1 /broad/dawnccle/sequencing/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/merged_fastqs/T47D-rep1_R1_bc_extracted.fastq.gz \
    -l /broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv \
    -o /broad/dawnccle/processed_data/missplicing_debug

"""

import os
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
    
    total_reads = 0
    aligned_reads = 0
    bc_not_found = 0
    included_reads_pre_chimera = 0
    chimera_reads = 0
    skipped_reads = 0
    
    unspliced_counter = 0
    minor_splice_counter = 0
    perfect_middle_exon = 0
    everything_passed_with_middle_exon = 0
    logging.info("Parsing fastq files.")

    # Create a new dictionary to store the counts with splice information.
    element_cb_umi_dict = {}
    chimera_reads_dict = {}
    # Parse the fastq. We only need to look at R1 since we already have the barcode and UMI.
    unzipped_file1 = gzip.open(fq1_file, "rt")
    for fq1 in SeqIO.parse(unzipped_file1, "fastq"):
        header = fq1.id
        read_name = header.split()[0]
        id, cb, umi = read_name.split("_")

        total_reads += 1
        # if total_reads % 10000 == 0:
        #     logging.info("Already processed reads: {}.".format(total_reads))

        # if total_reads >= 100000:
        #     break

        # Extract read and aligned region.
        guide_bc = cb
        cb_umi = cb + "_" + umi
        library_seq = str(fq1.seq)

        # We don't need to check anymore because the cell barcode is from the header whitelist.
        if guide_bc not in bc_lib_dict:
            bc_not_found += 1
            continue
        
        element_id = bc_lib_dict[guide_bc]
        if element_id != "ENSG00000170832.13;USP32;chr17-60207020-60207132-60205446-60205658-60208058-60208210":
            continue
        
        # Find the start and end of the "upstream exon" and "downstream exon" in the library sequence.
        # We will use this to determine if the read is spliced or not.
        start_upstream, end_upstream = find_substring_with_mismatches(
            library_seq, UPSTREAM_EXON_LAST_20, 1
        )

        # If the upstream exon is not found, then we can't do anything.
        # We are going to skip this, there is 4 out of 10,000 reads that don't have this.
        if start_upstream is None or end_upstream is None:
            continue

        # Trim the library sequence to exclude the upstream exon now that we have it mapped.
        library_seq_upstream_exon_trimmed = library_seq[end_upstream:]

        ########################################
        # Check out unspliced reads.
        ########################################
        # Get the first 20 or whatever is left.
        if len(library_seq_upstream_exon_trimmed) < 20:
            library_seq_first_20 = library_seq_upstream_exon_trimmed
        else:
            library_seq_first_20 = library_seq_upstream_exon_trimmed[:20]
        
        # Try to map to the unspliced reference.
        unsplice_start, _ = find_substring_with_mismatches(
            UNSPLICED_INTRON, library_seq_first_20, 1
        )
        if unsplice_start is not None:
            aligned_reads += 1
            # The NA at the end specifies that there is no inclusion length.
            temp_ID = f"{element_id}__UNSPLICED__{unsplice_start}__NA"
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
        # First we try to find the downstream exon.
        start_downstream, end_downstream = find_substring_with_mismatches(
            library_seq_upstream_exon_trimmed, DOWNSTREAM_EXON_FIRST_20, 0
        )
        
        # First assign the coordinate of the start of the downstream exon to 0 (0 means that the downstream exon splice site is what we expected, without any offset).
        # Basically, if the start_downsteam is found, then that means ref_downstream_exon_start is 0.
        ref_downstream_exon_start = -9999
        if start_downstream is not None and end_downstream is not None:
            ref_downstream_exon_start = 0
        
        # If the downstream exon is not found, we are going to loop through the "constant downstream exon" and see if we can find it. 
        # This is because maybe there is another splice site in the downstream exon.
        if start_downstream is None or end_downstream is None:
            ref_downstream_exon_found = False
            for i in range(len(DOWNSTREAM_INTRON) + len(DOWNSTREAM_EXON)):
                temp_substring = DOWNSTREAM_INTRON + DOWNSTREAM_EXON[i:i+20]
                start_downstream, end_downstream = find_substring_with_mismatches(
                    library_seq_upstream_exon_trimmed, temp_substring, 0
                )
                if start_downstream is not None and end_downstream is not None:
                    ref_downstream_exon_found = True
                    ref_downstream_exon_start = i - len(DOWNSTREAM_INTRON)
                    break
            # We skip if there is no downstream exon found.
            if not ref_downstream_exon_found:
                continue
        
        # If start_downstream is <=2, we keep it as a skipped read and just include the offset. 
        if start_downstream <= 2:
            aligned_reads += 1
            skipped_reads += 1
            temp_ID = f"{element_id}__SKIPPED__{ref_downstream_exon_start}__{start_downstream}"
            element_cb_umi_dict = add_cb_umi_to_dict(
                element_cb_umi_dict, temp_ID, cb_umi
            )
            continue    
            
        ########################################
        # Check out included sequences.
        ########################################
        # Get the reference sequence.
        ref_included_upstream_intron = upstream_intron_fasta[element_id]
        ref_included_downstream_exon = downstream_intron_fasta[element_id]
        ref_included_middle_exon = middle_exon_fasta[element_id]
        ref_included_sequence = ref_included_upstream_intron + ref_included_middle_exon + ref_included_downstream_exon + DOWNSTREAM_INTRON
        
        # Get the middle included sequence.
        middle_included_sequence = library_seq_upstream_exon_trimmed[:start_downstream]
        # Try to find this middle sequence in the reference sequence.
        # If <=20, we allow 0 mismatches.
        if len(middle_included_sequence) <= 20:
            start_middle_included_sequence, end_middle_included_sequence = find_substring_with_mismatches(
                ref_included_sequence, middle_included_sequence, 0
            )
        # If between 20 and 40, we allow 1 mismatch.
        elif len(middle_included_sequence) > 20 and len(middle_included_sequence) <= 40:
            start_middle_included_sequence, end_middle_included_sequence = find_substring_with_mismatches(
                ref_included_sequence, middle_included_sequence, 1
            )
        else:
            start_middle_included_sequence, end_middle_included_sequence = find_substring_with_mismatches(
                ref_included_sequence, middle_included_sequence, 2
            )
        # Now we can find the coordinate of the included sequence.
        aligned_reads += 1
        included_reads_pre_chimera += 1
        # If the middle included sequence is not found, we skip.
        if start_middle_included_sequence is None or end_middle_included_sequence is None:
            chimera_reads += 1
            continue
        
        start_guide_pos_adjusted = start_middle_included_sequence - len(ref_included_upstream_intron)
        end_guide_pos_adjusted = end_middle_included_sequence - len(ref_included_upstream_intron) - len(ref_included_middle_exon)
        
        temp_ID = f"{element_id}__INCLUDED__{start_guide_pos_adjusted}:{end_guide_pos_adjusted}:{ref_downstream_exon_start}__{len(middle_included_sequence)}"
        # if start_guide_pos_adjusted != 0 and end_guide_pos_adjusted != 0:
        print(f"{temp_ID},{middle_included_sequence},{library_seq}")
        element_cb_umi_dict = add_cb_umi_to_dict(element_cb_umi_dict, temp_ID, cb_umi)
        
        if start_guide_pos_adjusted == 0 and end_guide_pos_adjusted == 0:
            perfect_middle_exon += 1
        everything_passed_with_middle_exon += 1

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
    
    # Chimera reads dict.
    chimera_reads_df = pd.DataFrame.from_dict(chimera_reads_dict, orient="index")

    # Log all the read statistics and create CSV log entries
    stats = {
        "total_reads": total_reads,
        "total_aligned_reads": aligned_reads, 
        "bc_not_found": bc_not_found,
        "perc_bc_not_found": bc_not_found / float(total_reads),
        "skipped_reads": skipped_reads,
        "chimera_reads": chimera_reads,
        "perc_chimera_reads": chimera_reads / float(included_reads_pre_chimera),
        "unspliced_counter": unspliced_counter,
        "minor_splice_counter": minor_splice_counter,
        "perfect_middle_exon": perfect_middle_exon,
        "everything_passed_with_middle_exon": everything_passed_with_middle_exon
    }

    # Log all statistics and build stats_log array
    stats_log = []
    for stat_name, stat_value in stats.items():
        logging.info(f"{stat_name}: {stat_value}")
        stats_log.append(f"{stat_name},{stat_value}\n")

    return (element_cb_umi_df, umi_dedup_df, stats_log, chimera_reads_df)


def read_barcode_to_lib_table(table_path):
    d = {}
    f = open(table_path)
    lines = f.readlines()
    for line in lines:
        # Skip the header line.
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

    logging.info(
        "################# Starting on file: {} #################".format(
            Path(fq1_path)
        )
    )

    # Read barcode to element matching.
    library_bc_path = args.library_bc_path
    logging.info("Reading barcode from library table: {}".format(library_bc_path))
    barcode_lib_dict = read_barcode_to_lib_table(library_bc_path)

    # This is the processed dataframe that's important.
    identity_df, umi_dedup_df, stats_log, chimeric_reads_df = get_identity_from_fq(
        fq1_path,
        fq2_path,
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
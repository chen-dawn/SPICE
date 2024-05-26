'''
python /broad/thechenlab/Jing/find_misspliced.py \
    -f /broad/thechenlab/NextSeqOutput_Xerneas/230426_VL00297_175_AACM7VVM5/Data/Intensities/BaseCalls/K562_WT_1_S1_R1_bc_extracted.fastq.gz \
    -l /broad/thechenlab/Jing/splicing_analysis/20230511_library_cb_rev_comp.csv \
    -o /broad/thechenlab/Jing/splicing_output/230520_K562_misspliced_length_test
'''



import csv
import argparse
from Bio import SeqIO
import re
import pandas as pd
import logging
import sys
import os
from pathlib import Path
import statistics
import gzip


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
        "-f", "--fq1", type=str, required=True, help="fastq file",
    )

    parser.add_argument(
        "-l",
        "--library_bc_path",
        type=str,
        required=True,
        help="The location of the barcode/library matching dictionary.",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help="Directory where output files will be saved.",
        default="/broad/thechenlab/Jing/splicing_output/",
    )

    parser.add_argument(
        "-c",
        "--cell_barcode_path",
        type=str,
        help="The location of the cell barcode list as txt file.",
        default="/broad/thechenlab/Jing/splicing_output_230516/230619_K562_misspliced_seq/cb_splicing_insertions_results.txt",
    )

    args = parser.parse_args()
    return args


def find_substring_with_mismatches(string, substring, max_mismatches):
    substring_length = len(substring)
    max_index = len(string) - substring_length

    for i in range(max_index + 1):
        mismatches = 0

        for j in range(substring_length):
            if string[i+j] != substring[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break

        if mismatches <= max_mismatches:
            return i

    return -1


def make_exon_dict(table_path):
    exon_dict = {}
    table = csv.reader(open(table_path, 'r'), delimiter=',')
    next(table)
    for row in table:
        exon_dict[row[1]] = row[3]

    return exon_dict


def dict_mode(length_dict, seq_dict):
    new_length_dict = {}
    new_seq_dict = {}
    for key, values in length_dict.items():
        mode = statistics.mode(values)
        new_length_dict[key] = mode
        for seq in seq_dict[key]:
            if len(seq) == mode:
                new_seq_dict[key] = seq
                break

    return new_length_dict, new_seq_dict  


def read_cb_from_txt(file):
    cb_list = []
    with open(file) as f:
        for line in f:
            cb_list.append(line.strip())
    return cb_list
    

def get_misspliced(fq1_file, exon_dict, cb_list):

    up_insertion_count = {}
    down_insertion_count = {} 
    up_insertion_length = {}
    down_insertion_length = {}

    up_insertion_seq = {}

    upstream_exon = 'ATGATGTGGACCTGGAACAG'  # last 20 bp
    downstream_exon = 'GTGCGGCAGCTGGTGCCTCG' # first 20 bp

    input_file = gzip.open(fq1_file,"rt")

    total_reads = 0

    upstream_insertion = 0
    downstream_insertion = 0
    no_upstream_match = 0
    no_downstream_match = 0

    # match the name of the sequence to the library barcode in the dictionary
    for record in SeqIO.parse(input_file, "fastq"):
        total_reads += 1
        if total_reads % 10000 == 0:
            logging.info("Already processed reads: {}.".format(total_reads))
        
        header, seq = record.id, str(record.seq)
        read_name = header.split()[0]
        id, cb, umi = read_name.split("_")

        # check if the read match the cb list
        if cb not in cb_list:
            continue

        # if total_reads > 100000:
        #     break

        cb_umi = cb + "_" + umi

        library_exon = exon_dict[cb]
        # print(element_id)
        # print('seq:')
        # print(seq)
        # print('library exon:')
        # print(library_exon)

        up_exon_pos = find_substring_with_mismatches(seq, upstream_exon, 1)
        lib_exon_pos_1 = find_substring_with_mismatches(seq, library_exon[:20], 1)
        lib_exon_pos_2 = find_substring_with_mismatches(seq, library_exon[-20:], 1)
        down_exon_pos = find_substring_with_mismatches(seq, downstream_exon, 1)

        # up_exon_pos = seq.find(upstream_exon)
        # lib_exon_pos_1 = seq.find(library_exon[:20])
        # lib_exon_pos_2 = seq.find(library_exon[-20:])
        # down_exon_pos = seq.find(downstream_exon)


        if up_exon_pos == -1 or lib_exon_pos_1 == -1:
            no_upstream_match += 1

        if lib_exon_pos_2 == -1 or down_exon_pos == -1:
            no_downstream_match += 1


        if up_exon_pos != -1 and lib_exon_pos_1 != -1:
            up_dist = lib_exon_pos_1 - up_exon_pos - 20

            # if up_dist < 0:
            #     print(cb_umi)
            #     print('up_exon_pos:{}' .format(up_exon_pos))
            #     print('lib_exon_pos_1:{}' .format(lib_exon_pos_1))
            #     print('library_exon:{}'.format(library_exon))
            #     print(up_dist)
            #     print(seq)

            if up_dist > 1:
                upstream_insertion += 1
  
            if cb_umi not in up_insertion_length:
                up_insertion_length[cb_umi] = [up_dist]
            else:
                up_insertion_length[cb_umi].append(up_dist)
            
            if cb_umi not in up_insertion_count:
                up_insertion_count[cb_umi] = 1
            else:
                up_insertion_count[cb_umi] += 1

            if cb_umi not in up_insertion_seq:
                up_insertion_seq[cb_umi] = [seq[up_exon_pos+20:lib_exon_pos_1]]
            else:
                up_insertion_seq[cb_umi].append(seq[up_exon_pos+20:lib_exon_pos_1])
        
            up_insertion_length

    up_insertion_length_mode, up_insertion_seq_mode = dict_mode(up_insertion_length, up_insertion_seq)
        
    merged_df = pd.DataFrame(
        {  
            "upstream_insertion_count": pd.Series(up_insertion_count),
            "upstream_insertion_length": pd.Series(up_insertion_length_mode),
            "upstream_insertion_seq": pd.Series(up_insertion_seq_mode)
        }
    )

    logging.info("Total reads: {}".format(total_reads))

    return merged_df


if __name__ == "__main__":
    args = parse_args()
    fq1_path = args.fq1
    out_dir = args.output_dir
    cb_list = read_cb_from_txt(args.cell_barcode_path)

    # Create output dir if it doesn't exist.
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Set up logging.
    logging_filename = os.path.basename(fq1_path).replace("R1_bc_extracted.fastq.gz", "misspliced_analysis.log")
    logging_filename = os.path.join(out_dir, logging_filename)

    # Start logging.
    start_logger(logging_filename)
    logging.info("################# Starting on file: {}".format(Path(fq1_path)))
    
    # Read barcode to element matching.
    library_bc_path = args.library_bc_path
    logging.info("Reading barcode from library table: {}".format(library_bc_path))
    lib_exon_dict = make_exon_dict(library_bc_path)


    # This is the processed dataframe that's important.
    misspliced_df= get_misspliced(fq1_path, lib_exon_dict, cb_list)

    # save the df to a csv file
    csv_filename = os.path.basename(fq1_path).replace("R1_bc_extracted.fastq.gz", "misspliced_analysis.csv")
    csv_filename = os.path.join(out_dir, csv_filename)
    logging.info("Writing output table to: {}".format(csv_filename))
    misspliced_df.to_csv(csv_filename)
    


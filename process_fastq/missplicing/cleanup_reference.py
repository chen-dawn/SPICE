"""
Clean up the reference fasta file. 

Basically I realized that there are some sequeneces that lead to multimapping. Which is not good. 
So we are going to deduplicate the reference fasta file. We will do this by first sorting the middle exon from the shortest to the longest. 
Then, we will check if the middle exon is uniquely mapped. If it is, we will add it to the definitely unique dict.
If it's multimapped, we will keep the shortest sequence and then add the rest to the discard pile if there are 100% perfectly maapped.
This leads to a discard of ~1.7k from the reference, which is alright.

Number of definitely unique middle exons: 41137
Number of multi-mapped middle exons: 3511
Number of multi-mapped guides to discard: 1724
Number of not mapped middle exons: 0

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

guide_fasta_path = (
    "/broad/dawnccle/melange/data/guide_library/PRISM_47k_twist.fasta"
)
middle_exon_fasta_path = (
    "/broad/dawnccle/melange/data/guide_library/PRISM_47k_skipped_exon.fasta"
)
upstream_intron_fasta_path = (
    "/broad/dawnccle/melange/data/guide_library/PRISM_47k_upstream_intron.fasta"
)
downstream_intron_fasta_path = (
    "/broad/dawnccle/melange/data/guide_library/PRISM_47k_downstream_intron.fasta"
)


middle_exon_fasta = read_guide_fasta(middle_exon_fasta_path)
upstream_intron_fasta = read_guide_fasta(upstream_intron_fasta_path)
downstream_intron_fasta = read_guide_fasta(downstream_intron_fasta_path)
guide_fasta = read_guide_fasta(guide_fasta_path)


aln = mp.Aligner(
        fn_idx_in=guide_fasta_path, preset="sr", best_n=10,k=9,w=1,min_chain_score=20)

# Convert the middle exon fasta to a list of tuples.
middle_exon_list = list(middle_exon_fasta.items())
# Sort the middle exon list by the length of the sequence which is the 2nd column. from smallest to largest.
middle_exon_list.sort(key=lambda x: len(x[1]))

already_considered = set()
definitely_unique_dict = {}
multi_mapped_to_include_dict = {}
multi_mapped_to_discard_dict = {}
alt_ref_dict = {}
not_mapped_count = 0
for i, (middle_exon_id, middle_exon_seq) in enumerate(middle_exon_list):
    if middle_exon_id in already_considered:
        continue
    if i % 1000 == 0:
        print(f"Processing middle exon {i}/{len(middle_exon_list)}")
    num_mapped = 0
    # print("##################")
    for al in aln.map(middle_exon_seq):
        # print(al)
        mapped_guide_id = al.ctg
        # Get the guide id.
        num_mapped += 1
    
    if num_mapped == 1:
        # Sanity check that the al.ctg is the same.
        if mapped_guide_id != middle_exon_id:
            logging.error(f"Middle exon {middle_exon_id} mapped to {mapped_guide_id}")
            continue
        
        definitely_unique_dict[middle_exon_id] = definitely_unique_dict.get(middle_exon_id, 0) + 1
        # Add to the already considered set.
        already_considered.add(middle_exon_id)
        
    elif num_mapped > 1:
        # Because we are doing from shortest to longest, we can assume that the middle_exon_id is the best match.
        multi_mapped_to_include_dict[middle_exon_id] = multi_mapped_to_include_dict.get(middle_exon_id, 0) + 1
        already_considered.add(middle_exon_id)
        for al in aln.map(middle_exon_seq):
            mapped_guide_id = al.ctg
            if mapped_guide_id != middle_exon_id:
                # We are only interested in the multi-mapped guides that are perfect matches at the most length.
                if al.mlen == len(middle_exon_seq) and al.NM == 0:
                    multi_mapped_to_discard_dict[mapped_guide_id] = multi_mapped_to_discard_dict.get(mapped_guide_id, 0) + 1
                    alt_ref_dict[mapped_guide_id] = middle_exon_id
                    already_considered.add(mapped_guide_id)
            
    else:
        print("NOT MAPPED?")
        print(middle_exon_id)
        not_mapped_count += 1
        print(not_mapped_count)

print(f"Number of definitely unique middle exons: {len(definitely_unique_dict)}")
print(f"Number of multi-mapped middle exons: {len(multi_mapped_to_include_dict)}")
print(f"Number of multi-mapped guides to discard: {len(multi_mapped_to_discard_dict)}")
print(f"Number of not mapped middle exons: {not_mapped_count}")

# Write to file.
f = open("/broad/dawnccle/melange/data/guide_library_cleaned/ref_test_definitely_unique_middle_exons.tsv", "w")
for k, v in definitely_unique_dict.items():
    f.write(f"{k}\t{v}\n")
    
f = open("/broad/dawnccle/melange/data/guide_library_cleaned/ref_test_multi_mapped_middle_exons_to_include.tsv", "w")
for k, v in multi_mapped_to_include_dict.items():
    f.write(f"{k}\t{v}\n")
    
f = open("/broad/dawnccle/melange/data/guide_library_cleaned/ref_test_multi_mapped_guides_to_discard.tsv", "w")
for k, v in multi_mapped_to_discard_dict.items():
    f.write(f"{k}\t{v}\n")
    
f = open("/broad/dawnccle/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv", "w")
f.write("alt\tref\n")
for k, v in alt_ref_dict.items():
    f.write(f"{k}\t{v}\n")
# How the refernece is generated

## Map the duplicate regions to a reference region
First I ran the script `/Volumes/broad_dawnccle/melange/process_fastq/missplicing/cleanup_reference.py` which finds the sequences that are definitely unique. And for the sequences that are not unique, we keep the shortest one and discard the longer ones. This generates a list of mapping tables:

```
/broad/dawnccle/melange/data/guide_library_cleaned/ref_test_definitely_unique_middle_exons.tsv
/broad/dawnccle/melange/data/guide_library_cleaned/ref_test_multi_mapped_middle_exons_to_include.tsv
/broad/dawnccle/melange/data/guide_library_cleaned/ref_test_multi_mapped_guides_to_discard.tsv
```

And also a dictionary that has the alt to ref exons.
```
/broad/dawnccle/melange/data/guide_library_cleaned/ref_test_alt_ref_dict.tsv
```

## Filter the reference to remove duplicate sequences
I also run the script `/Volumes/broad_dawnccle/melange/process_fastq/missplicing/clean_up_reference.R` that takes in the reference and updates the reference to include the alt exons. This generates a new reference file:

```
/Volumes/broad_dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv
/Volumes/broad_dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_filtered.csv
```

These files are then used to generate the reference for the guide library and are filtered. 

## Generate the reference fasta files
Finally, I run the script `/Volumes/broad_dawnccle/melange/data/guide_library_cleaned/MakeRNAElementDict_47klib_updated.py` which generates the reference fasta files.
The files are:
```
    guide_fasta_out_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_guides_filtered.fasta"
    skipped_exon_fasta_out_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_skipped_exon_filtered.fasta"
    upstream_intron_fasta_out_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_upstream_intron_filtered.fasta"
    downstream_intron_fasta_out_path = "/broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_downstream_intron_filtered.fasta"
    library_reference_dir = /broad/dawnccle/melange/data/guide_library_cleaned/WEAK_47k_reference_no_adapter_filtered.fasta"
```

# How the fastq files are processed

## Use UMI_tools to extract the element barcode from R2
First I run umi tools using the barcode whitelist like so to extract the barcode from R2. 
This is the file `/Volumes/broad_dawnccle/melange/process_fastq/V1_original_pipeline/run_umi_tools_demux.sh`. 

The driver function is:
```
# Run like:
# while read p; do
#     echo $p
#     qsub -v FILENAME=$p /Volumes/broad_dawnccle/melange/process_fastq/V1_original_pipeline/run_umi_tools_demux.sh
# done < merged_R1_filelist.txt

cd /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/merged_fastqs

echo $FILENAME
BASENAME=${FILENAME%_merged_R1_001.fastq.gz}
FQ1=${BASENAME}_merged_R1_001.fastq.gz
FQ2=${BASENAME}_merged_R2_001.fastq.gz

umi_tools extract --stdin $FQ1 --read2-in=${FQ2}  \
    --bc-pattern2="(?P<umi_1>.{10})(?P<discard_2>TTGCTAGGACCGGCCTTAAAGC){s<=2}(?P<cell_1>.{14})(?P<discard_3>CTACTTA){s<=1}.*" \
    --extract-method=regex \
    --whitelist /broad/thechenlab/Dawn/splicing/full_library_unique_exons_only/twist_47k_barcode_whitelist_reverse_complement.txt  \
    -L ${BASENAME}.log --stdout ${BASENAME}_R1_bc_extracted.fastq.gz --read2-out ${BASENAME}_R2_bc_extracted.fastq.gz \
    --filtered-out ${BASENAME}_R1_bc_extracted_failed.fastq.gz --filtered-out2 ${BASENAME}_R2_bc_extracted_failed.fastq.gz
```

## Get the splicing information from R1 and R2
I run the following scripts to call splicing information from each read.

```
# Run for all samples in the big screen. 
/Volumes/broad_dawnccle/melange/process_fastq/missplicing/run_barcode_match_to_element.sh
# K562 WT and K700E only, collected separately. 
/Volumes/broad_dawnccle/melange/process_fastq/missplicing/run_barcode_match_to_element_K562.sh
```

The driver function that is in these script is:
```
python /broad/dawnccle/melange/process_fastq/missplicing/MatchBarcodeToElementRNA_FINAL_VERSION_V5_240715.py \
    -1 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/merged_fastqs/${BASENAME}_R1_bc_extracted.fastq.gz \
    -l /broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv \
    -o /broad/dawnccle/processed_data/missplicing_test_v6
```

Overall, we get 3 output files from each sample. Just as an example:
```
# The element barcode umi count matrix for every UMI. this file is huge.
T47D-rep1_element_cb_umi_count_fine_grained_idx.csv
# The stats log file for the file, contains information like how many reads mapped and what the alignment states are.
T47D-rep1_stats_log_fine_grained_idx.txt
# The element, barcode, umi count matrix for every splicing event. This is the file that is used for downstream analysis.
T47D-rep1_umi_dedup_fine_grained_idx.csv
```

*The most current latest version is always saved in the `/broad/dawnccle/processed_data/latest` folder.*

## Merge all the files together
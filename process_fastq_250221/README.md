# Process and normalize FASTQ files

## 1. Convert raw FASTQ files to counts table

### Use UMI-tools to extract the element barcode from R2
Our sequencing is such that R1 contains the information for the reporter sequence and splice junctions, and R2 contains the element barcode. First I run umi tools using the barcode whitelist like so to extract the barcode from R2. This only allows for perfect string matching with no error tolerance, but the Illumina sequencing error rate is low enough such that it doesn't lead to a major decrease in mapping rates.

The file to run is `/Volumes/broad_dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/run_umi_tools_demux.sh`. 

The driver function within the file is:
```
# Run like:
# merged_R1_filelist.txt is a list of all the R1 filenames.
# while read p; do
#     echo $p
#     qsub -v FILENAME=$p /Volumes/broad_dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/run_umi_tools_demux.sh
# done < merged_R1_filelist.txt

cd /broad/dawnccle/sequencing/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/merged_fastqs

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

### Get the splicing information from R1 and R2
I run the following scripts to call splicing information from each read.

```
# Run for all samples in the big screen. 
/broad/dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/run_barcode_match_to_element_nova230516.sh
# K562 WT and K700E only, collected separately. 
/broad/dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/run_barcode_match_to_element_K562.sh
# Run for Nova 240826 and 241106.
/broad/dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/run_barcode_match_to_element_nova240826.sh
/broad/dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/run_barcode_match_to_element_nova241106.sh
```

The driver function that is in these script is:
```
python /broad/dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/MatchBarcodeToElementRNA_20250220.py \
    -1 ${indir}/${BASENAME}_R1_bc_extracted.fastq.gz \
    -l /broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv \
    -o $outdir
```

We are mapping it to the new, cleaned reference library `/broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv`. This library has 44648 sequences (instead of 46372 in the full). This is because I found that there were some sequences that had the "same" skipped exon, just that the annotation positions were slightly different. So previously, I have filtered out the sequences that had the same skipped exons. 

In the script, we are using the following reference files:
```
    /broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_guides_filtered.fasta
    /broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_skipped_exon_filtered.fasta
    /broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_upstream_intron_filtered.fasta
    /broad/dawnccle/melange/data/guide_library_cleaned/PRISM_47k_downstream_intron_filtered.fasta
```

Overall, we get 3 output files from each sample. Just as an example:
```
# The element barcode umi count matrix for every UMI. this file is huge.
T47D-rep1_element_cb_umi_count_fine_grained_idx.csv

# The stats log file for the file, contains information like how many reads mapped and what the alignment states are.
T47D-rep1_stats_log_fine_grained_idx.txt

# The element, barcode, umi count matrix for every splicing event. 
# This is the file that is used for downstream analysis.
T47D-rep1_umi_dedup_fine_grained_idx.csv
```
** The most recent set of files are saved in `/broad/dawnccle/processed_data/reprocess_250221/` **


## 2. Merge and normalize the counts table

### Refer to the README in the `02_merge_and_normalize_counts` folder for more details.

## 3. Convert to PSI

### Refer to the README in the `03_convert_to_PSI` folder for more details.



















# Statistics of differential splicing using the rMATS-STAT pipeline

Here we use the [rMATS-STAT](https://github.com/Xinglab/rMATS-STAT) pipeline to calculate the differential splicing statistics. I made some changes to this script to remove most of the print statements. Some of the files are named with "Pairiadise" but that's incorrect, just ignore - we are really running the standalone rMATS-STAT package and I just didn't change the name of the scripts. 

## Comparing one cell type against all other cell types
The scripts we are running are:
```
/Volumes/broad_dawnccle/melange/process_fastq/final_pipeline/run_pairadise_pair_indiv_PSI_one_to_all.sh
/Volumes/broad_dawnccle/melange/process_fastq/final_pipeline/run_pairadise_pair_indiv_3ss_one_to_all.sh
```

I am running it as follows:
```
# Read the cell types into an array
readarray -t celltypes </broad/dawnccle/melange/process_fastq/missplicing/unique_celltypes_names_only.txt

# Loop over each pair of cell types
for ((i = 0; i < ${#celltypes[@]}; i++)); do
    celltype1=${celltypes[i]}
    echo "Running for cell type: $celltype1"
    qsub -v celltype1=$celltype1 /broad/dawnccle/melange/process_fastq/final_pipeline/run_pairadise_pair_indiv_PSI_one_to_all.sh
    qsub -v celltype1=$celltype1 /broad/dawnccle/melange/process_fastq/final_pipeline/run_pairadise_pair_indiv_3ss_one_to_all.sh
done
```

## Comparing one tissue type against another tissue type
Same as above, we just run one tissue type against another tissue type. 
```
# Read the cell types into an array
readarray -t celltypes </broad/dawnccle/melange/data/unique_lineages.txt

# Create a set to track processed pairs
declare -A processed_pairs
# Reset the processed_pairs array
processed_pairs=()

# Loop over each pair of cell types
for ((i = 0; i < ${#celltypes[@]}; i++)); do
    for ((j = i + 1; j < ${#celltypes[@]}; j++)); do
        celltype1=${celltypes[i]}
        celltype2=${celltypes[j]}

        # Skip if celltype1 is the same as celltype2
        if [[ "$celltype1" == "$celltype2" ]]; then
            continue
        fi

        # Create a sorted key to ensure pairs are unique
        pair_key=$(echo "$celltype1 $celltype2" | tr " " "\n" | sort | tr "\n" "_")

        # Skip if the pair has already been processed
        if [[ -n "${processed_pairs[$pair_key]}" ]]; then
            continue
        fi

        # Mark the pair as processed
        processed_pairs[$pair_key]=1

        echo "Running for pair: $celltype1 and $celltype2"
        qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq/final_pipeline/run_pairadise_pair_tissuetype_PSI.sh
        qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq/final_pipeline/run_pairadise_pair_tissuetype_3ss.sh
    done
done
```

## Combine the files for all the individual comparisons

Run the following for each of the comparison folders. This only filters for values that have <0.01 FDR.
```
output_file="rmats_one_vs_all_combined_output_PSI.tsv"
# output_file="rmats_one_vs_all_combined_output_3ss.tsv"
# output_file="rmats_tissue_type_combined_output_PSI.tsv"
# output_file="rmats_tissue_type_combined_output_3ss.tsv"
first_file=true
for folder in */; do
    folder_name="${folder%/}"
    for file in "$folder"/*FDR.txt; do
        echo -e "\n####################"
        echo "$file"
        echo -e "####################\n"

        if $first_file; then
            awk 'NR==1 {print $0 "\tFolder"} NR>1 && $9 < 0.01 {print $0 "\t\"'$folder_name'\""}' "$file" >>"$output_file"
            first_file=false
        else
            awk 'NR>1 && $9 < 0.01 {print $0 "\t\"'$folder_name'\""}' "$file" >>"$output_file"
        fi
    done
done
```
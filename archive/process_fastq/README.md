# Process FASTQ files to generate count matrix

**The most current latest version for all processed files is always saved in the `/broad/dawnccle/processed_data/latest` folder.**

## Use UMI-tools to extract the element barcode from R2
Our sequencing is such that R1 contains the information for the reporter sequence and splice junctions, and R2 contains the element barcode. First I run umi tools using the barcode whitelist like so to extract the barcode from R2. This only allows for perfect string matching with no error tolerance, but the Illumina sequencing error rate is low enough such that it doesn't lead to a major decrease in mapping rates.

The file to run is `/Volumes/broad_dawnccle/melange/process_fastq/V1_original_pipeline/run_umi_tools_demux.sh`. 

The driver function within the file is:
```
# Run like:
# merged_R1_filelist.txt is a list of all the R1 filenames.
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
    -l /broad/dawnccle/melange/data/guide_library/20230130_twist_library_v3_ID_barcode_ROUT.csv \
    -o /broad/dawnccle/processed_data/missplicing_test_v6
```

Here, we still map it to the old, full barcode reference `/broad/dawnccle/melange/data/guide_library/20230130_twist_library_v3_ID_barcode_ROUT.csv` even though I removed some of the sequences from the sequence reference, because in the subsequent merging step (using `01_merge_counts_all_samples_workstation.R`), I will also merge the "exon skipped" barocdes to the consensus reference. Since for the "included" case, the aligner will already just align it to the correct sequence, not mapping back the exon skipped sequences to the reference would result in undercounting and inflated inclusion rates. 

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

I also pretty clean the data using the following script (it's not really necessary, but I happen to be using this "cleaned" file for downstream processing):
```
# To pretty clean the files.
for file in *umi_dedup_fine_grained_idx.csv; do
    echo $file
    awk -F',' 'NR==1 {print "index\tmode\toffset\tcount"} 
               NR>1 {
                   split($1, id_parts, "_")
                   name=""
                   for (i=1; i<=length(id_parts)-2; i++) {
                       name = (i == 1) ? id_parts[i] : name "_" id_parts[i]
                   }
                   mode = id_parts[length(id_parts)-1]
                   offset = id_parts[length(id_parts)]
                   print name"\t"mode"\t"offset"\t"$2
               }' $file >${file%.csv}_formatted.tsv
done
```

**The most current latest version is always saved in the `/broad/dawnccle/processed_data/latest` folder.**
Here are the paths of where the latest files are stored:

```
# All samples (same as file in V5/V6 processed files)
/broad/dawnccle/processed_data/latest/raw_47celltype

# K562 WT vs K700E (same as files in V6 processed files, accidentially deleted one of the V5 files)
/broad/dawnccle/processed_data/latest/raw_K700E
```

## Merge all the files together
Run the script `/broad/dawnccle/melange/process_fastq/final_pipeline/01_merge_counts_all_samples_workstation.R`. This merges all the individual files and saves them to the directory.

The main files that should be used are:
```
/broad/dawnccle/processed_data/latest/umi_count_merged_to_ref_normalized.csv
/broad/dawnccle/processed_data/latest/K700E_umi_count_merged_to_ref_normalized.csv
```

This also makes the files that have some minor filtering which are:
```
# This one is for PSI calculations, where only the reads that are perfectly skipped or included (as according to the reference file) are kept.
/broad/dawnccle/processed_data/latest/all_sample_reps_PSI.csv
# This one is for alternative 5'ss calculations.
/broad/dawnccle/processed_data/latest/all_sample_reps_3ss.csv
```

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
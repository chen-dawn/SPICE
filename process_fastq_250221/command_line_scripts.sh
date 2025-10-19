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

######################
# This is for celltype pairs
######################
# Read the cell types into an array
readarray -t celltypes </broad/dawnccle/melange/process_fastq/missplicing/unique_celltypes.txt

# Create a set to track processed pairs
declare -A processed_pairs

# Specify the output file for the qsub commands
output_file="qsub_commands_PSI.txt"

# Ensure the output file is empty at the start
>"$output_file"

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
        # Write the qsub command to the output file
        echo "qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq/missplicing/run_pairadise_pair_celltype_PSI.sh" >>"$output_file"
    done
done

######################
# This is for tissue type pairs
######################
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
        qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq/missplicing/run_pairadise_pair_tissuetype.sh
        # sh /broad/dawnccle/melange/process_fastq/missplicing/run_pairadise_pair_tissuetype.sh $celltype1 $celltype2
    done
done

######################
# This is for tissue type pairs PSI
######################
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
        # sh /broad/dawnccle/melange/process_fastq/missplicing/run_pairadise_pair_tissuetype.sh $celltype1 $celltype2
    done
done

qstata | awk '{print $1}' | while read -r job_id; do
    echo "Deleting job: $job_id"
    qdel "$job_id"
done

# Merge all output files.
output_file="missplicing_indiv_combined_output.tsv"
first_file=true

for folder in */; do
    folder_name="${folder%/}"
    for file in "$folder"/*FDR.txt; do
        echo -e "\n####################"
        echo "$file"
        echo -e "####################\n"

        if $first_file; then
            awk 'NR==1 {print $0 "\tFolder"} NR>1 && $9 < 0.01 {print $0 "\t\"'$folder_name'\""}' "$file" | column -t >>"$output_file"
            first_file=false
        else
            awk 'NR>1 && $9 < 0.01 {print $0 "\t\"'$folder_name'\""}' "$file" | column -t >>"$output_file"
        fi
    done
done

output_file="PSI_combined_output.tsv"
first_file=true

for folder in */; do
    folder_name="${folder%/}"
    # if folder contains "blood" continue
    if [[ "$folder_name" == *"blood"* ]]; then
        continue
    fi
    for file in "$folder"/*FDR.txt; do
        echo -e "\n####################"
        echo "$file"
        echo -e "####################\n"

        awk 'NR==1 || $9 < 0.01 {print $0 "\t\"'$folder_name'\""}' "$file" | column -t

    done
done

# Merge all output files.
output_file="NovaseqPSI_indiv_combined_output.tsv"
first_file=true

for folder in */; do
    folder_name="${folder%/}"
    for file in "$folder"/*FDR.txt; do
        echo -e "\n####################"
        echo "$file"
        echo -e "####################\n"

        if $first_file; then
            awk 'NR==1 {print $0 "\tFolder"} NR>1 && $9 < 0.01 {print $0 "\t\"'$folder_name'\""}' "$file" | column -t >>"$output_file"
            first_file=false
        else
            awk 'NR>1 && $9 < 0.01 {print $0 "\t\"'$folder_name'\""}' "$file" | column -t >>"$output_file"
        fi
    done
done

######################
# This is for celltype each one only
######################
# Read the cell types into an array
readarray -t celltypes </broad/dawnccle/melange/process_fastq/missplicing/unique_celltypes_names_only.txt

# Specify the output file for the qsub commands
output_file="qsub_commands_PSI_indiv.txt"

# Ensure the output file is empty at the start
>"$output_file"

# Loop over each pair of cell types
for ((i = 0; i < ${#celltypes[@]}; i++)); do
    celltype1=${celltypes[i]}
    # Write the qsub command to the output file
    echo "qsub -v celltype1=$celltype1 /broad/dawnccle/melange/process_fastq/missplicing/run_pairadise_pair_celltype_PSI_one_to_all.sh" >>"$output_file"
done

######################
# This is for celltype each one only for 5'SS
######################
# Read the cell types into an array
readarray -t celltypes </broad/dawnccle/melange/process_fastq/missplicing/unique_celltypes_names_only.txt

# Specify the output file for the qsub commands
output_file="qsub_commands_indiv.txt"

# Ensure the output file is empty at the start
>"$output_file"

# Loop over each pair of cell types
for ((i = 0; i < ${#celltypes[@]}; i++)); do
    celltype1=${celltypes[i]}
    # Write the qsub command to the output file
    # echo "qsub -v celltype1=$celltype1 /broad/dawnccle/melange/process_fastq/missplicing/run_pairadise_pair_celltype_one_to_all.sh" >>"$output_file"
    qsub -v celltype1=$celltype1 /broad/dawnccle/melange/process_fastq/final_pipeline/run_pairadise_pair_indiv_3ss_one_to_all.sh
    qsub -v celltype1=$celltype1 /broad/dawnccle/melange/process_fastq/final_pipeline/run_pairadise_pair_indiv_PSI_one_to_all.sh
done






output_file="PSI_combined_output_indiv.tsv"
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

output_file="5ss_combined_output_indiv.tsv"
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


for FILENAME in *rep*_R1_001.fastq.gz; do
    echo $FILENAME
    BASENAME=${FILENAME%_R1_001.fastq.gz}
    FQ1=${BASENAME}_R1_001.fastq.gz
    FQ2=${BASENAME}_R2_001.fastq.gz

    umi_tools extract --stdin $FQ1 --read2-in=${FQ2} \
        --bc-pattern2="(?P<umi_1>.{10})(?P<discard_2>TTGCTAGGACCGGCCTTAAAGC){s<=2}(?P<cell_1>.{14})(?P<discard_3>CTACTTA){s<=1}.*" \
        --extract-method=regex \
        --whitelist /broad/thechenlab/Dawn/splicing/full_library_unique_exons_only/twist_47k_barcode_whitelist_reverse_complement.txt \
        -L ${BASENAME}.log --stdout ${BASENAME}_R1_bc_extracted.fastq.gz --read2-out ${BASENAME}_R2_bc_extracted.fastq.gz \
        --filtered-out ${BASENAME}_R1_bc_extracted_failed.fastq.gz --filtered-out2 ${BASENAME}_R2_bc_extracted_failed.fastq.gz

    # filebase=${filename%_R1_bc_extracted.fastq.gz}

    python /broad/thechenlab/Dawn/splicing/library_sequencing_47k/MatchBarcodeToElementRNA_umi_tools_extracted_47k_R1_VCP_230524.py \
        -1 ${BASENAME}_R1_bc_extracted.fastq.gz \
        -2 ${BASENAME}_R2_bc_extracted.fastq.gz \
        -l /broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv \
        -r /broad/thechenlab/Dawn/splicing/library_sequencing_47k/RNA_ref/WEAK_47k_reference.fasta \
        -o /broad/thechenlab/Dawn/NextSeq/240706_VL00297_10_AAFW7M2M5/Data/Intensities/BaseCalls/splicing_out
done

# New script, final script. 
for FILENAME in T47D*rep*_R1_001.fastq.gz; do
    echo $FILENAME
    BASENAME=${FILENAME%_R1_001.fastq.gz}
    FQ1=${BASENAME}_R1_001.fastq.gz
    FQ2=${BASENAME}_R2_001.fastq.gz

    python /broad/dawnccle/melange/process_fastq/missplicing/MatchBarcodeToElementRNA_FINAL_VERSION_V5_240715.py \
        -1 ${BASENAME}_R1_bc_extracted.fastq.gz \
        -l /broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv \
        -o /broad/thechenlab/Dawn/NextSeq/240706_VL00297_10_AAFW7M2M5/Data/Intensities/BaseCalls/splicing_out_new
done

# Old script that lets you output chimeric reads now. 
for FILENAME in *rep*_R1_001.fastq.gz; do
    echo $FILENAME
    BASENAME=${FILENAME%_R1_001.fastq.gz}
    FQ1=${BASENAME}_R1_001.fastq.gz
    FQ2=${BASENAME}_R2_001.fastq.gz

    python /broad/dawnccle/melange/process_fastq/missplicing/MatchBarcodeToElementRNA_umi_tools_extracted_Novaseq230524_missplicing.py \
        -1 ${BASENAME}_R1_bc_extracted.fastq.gz \
        -l /broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv \
        -o /broad/thechenlab/Dawn/NextSeq/240706_VL00297_10_AAFW7M2M5/Data/Intensities/BaseCalls/splicing_out_new
done

######################
# This is for celltype pairs Novaseq 240826
######################
# Read the cell types into an array
readarray -t celltypes </broad/dawnccle/melange/process_fastq/novaseq240826_pipeline/nova240826_celllines.txt

# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=()

# Specify the output file for the qsub commands
output_file="qsub_commands_PSI.txt"

# Ensure the output file is empty at the start
>"$output_file"

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
        # Write the qsub command to the output file
        echo "qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq/novaseq240826_pipeline/run_pairadise_pair_celltype_PSI.sh" >>"$output_file"
    done
done

# Merge the sequences for FDR
output_file="Nova240826_PSI_combined_output_indiv.tsv"
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

output_file="Nova240826_3ss_combined_output_indiv.tsv"
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

######################
# This is for celltype pairs Novaseq 241106
######################
# Read the cell types into an array
readarray -t celltypes </broad/dawnccle/melange/process_fastq/novaseq240826_pipeline/nova241106_RBP_celllines.txt

# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=()

# Specify the output file for the qsub commands
output_file="qsub_commands_PSI.txt"

# Ensure the output file is empty at the start
>"$output_file"

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
        # Write the qsub command to the output file
        echo "qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq/novaseq241106_pipeline/run_pairadise_pair_celltype_PSI.sh" >>"$output_file"
    done
done
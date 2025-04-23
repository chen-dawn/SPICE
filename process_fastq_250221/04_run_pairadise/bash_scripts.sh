
######################
# This is for OEx pairs
######################
# Compare the dox vs no dox only
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_PSI_OEx.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Loop over RBP identifiers (1 to 12)
for rbp in {1..12}; do
    # Define the patterns for dox and no_dox for the current RBP
    no_dox_celltype="splicelib_hek_no_dox_rbp${rbp}"
    dox_celltype="splicelib_hek_dox_rbp${rbp}"

        echo "Running for RBP${rbp}: $dox_celltype and $no_dox_celltype"
        # Write the qsub command to the output file
        echo "qsub -v celltype1=$dox_celltype,celltype2=$no_dox_celltype /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_OEx_PSI.sh" >>"$output_file"
done

######################
# This is for MUT pairs
######################
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_PSI_MUT.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Get all cell types in the MUT directory
celltypes=($(ls /broad/dawnccle/processed_data/reprocess_250221/count_normalized_v4/MUT/*.tsv | xargs -n 1 basename | sed -E 's/-rep[0-9]+.*//' | sort -u))
# Create a set to track processed pairs
declare -A processed_pairs
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
        # Write the qsub command to the output file
        echo "qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_MUT_PSI.sh" >>"$output_file"
    done
done


######################
# This is for WT pairs transcriptomic group
######################
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_PSI_WT_transcriptomic_group.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Get all cluster numbers
file_path=/broad/dawnccle/melange/data/transcriptomic_groups.txt
# Cluster number is just the first column of the file. Just get the first column.
cluster_numbers=($(cut -d',' -f1 "$file_path"))
# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=()

#  Loop over each cell type
for cluster_number in "${cluster_numbers[@]}"; do
    echo "Running for cluster number: $cluster_number"
    # Write the qsub command to the output file
    echo "qsub -v cluster_number=$cluster_number /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_WT_PSI_transcriptomic_group.sh" >>"$output_file"
done

######################
# This is for WT pairs one to all
######################
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_PSI_WT_one_to_all.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Get all cell types in the WT directory
celltypes=($(ls /broad/dawnccle/processed_data/reprocess_250221/count_normalized_v4/WT/*.tsv | xargs -n 1 basename | sed -E 's/-rep[0-9]+.*//' | sort -u))
# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=()

#  Loop over each cell type
for celltype in "${celltypes[@]}"; do
    echo "Running for celltype: $celltype"
    # Write the qsub command to the output file
    echo "qsub -v celltype1=$celltype,celltype2=all /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_WT_PSI.sh" >>"$output_file"
done


######################
# Merge the sequences for FDR for WT
######################
output_file="WT_PSI_combined_output_indiv.tsv"
first_file=true

for folder in WT_*all/; do
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
# Merge the sequences for FDR for MUT
######################
output_file="MUT_PSI_combined_output_indiv.tsv"
first_file=true

for folder in MUT_*/; do
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

# Merge the sequences for FDR for OEx
output_file="OEx_PSI_combined_output_indiv.tsv"
first_file=true

for folder in OEx_*/; do
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


# Merge the sequences for FDR for WT_transcriptomic_group
output_file="WT_transcriptomic_group_PSI_combined_output_indiv.tsv"
first_file=true

for folder in WT_transcriptomic_group_*/; do
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
    
# Merge the seqeunces for FDR for RBPs
output_file="RBP_PSI_combined_output_indiv.tsv"
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
# This is for WT pairs altSS one to all
######################
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_altSS_WT_one_to_all.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Get all cell types in the WT directory
celltypes=($(ls /broad/dawnccle/processed_data/reprocess_250221/count_normalized_v4/WT/*.tsv | xargs -n 1 basename | sed -E 's/-rep[0-9]+.*//' | sort -u))
# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=()

#  Loop over each cell type
for celltype in "${celltypes[@]}"; do
    echo "Running for celltype: $celltype"
    # Write the qsub command to the output file
    echo "qsub -v celltype1=$celltype,celltype2=all /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_WT_altSS.sh" >>"$output_file"
done

######################
# This is for MUT pairs altSS one to all
######################
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_altSS_MUT_one_to_all.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Get all cell types in the MUT directory
celltypes=($(ls /broad/dawnccle/processed_data/reprocess_250221/count_normalized_v4/MUT/*.tsv | xargs -n 1 basename | sed -E 's/-rep[0-9]+.*//' | sort -u))
# Create a set to track processed pairs
declare -A processed_pairs
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
        # Write the qsub command to the output file
        echo "qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_MUT_altSS.sh" >>"$output_file"
    done
done

######################
# This is for WT pairs transcriptomic group
######################
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_altSS_WT_transcriptomic_group.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Get all cluster numbers
file_path=/broad/dawnccle/melange/data/transcriptomic_groups.txt
# Cluster number is just the first column of the file. Just get the first column.
cluster_numbers=($(cut -d',' -f1 "$file_path"))
# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=()

#  Loop over each cell type
for cluster_number in "${cluster_numbers[@]}"; do
    echo "Running for cluster number: $cluster_number"
    # Write the qsub command to the output file
    echo "qsub -v cluster_number=$cluster_number /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_WT_altSS_transcriptomic_group.sh" >>"$output_file"
done

# Merge the sequences for FDR for WT.
output_file="WT_altSS_combined_output_indiv.tsv"
first_file=true

for folder in WT_*all/; do
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

# Merge the sequences for FDR for MUT.
output_file="MUT_altSS_combined_output_indiv.tsv"
first_file=true

for folder in MUT_*/; do
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

# Merge the sequences for FDR for transcriptomic group.
output_file="WT_transcriptomic_group_altSS_combined_output_indiv.tsv"
first_file=true

for folder in WT_transcriptomic_group_*/; do
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



######################
# This is for celltype pairs Novaseq 241106
######################
# Compare the dox vs no dox only
# Specify the output file for the qsub commands
output_file="qsub_commands_PSI_RBP.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Loop over RBP identifiers (1 to 12)
for rbp in {1..12}; do
    # Define the patterns for dox and no_dox for the current RBP
    no_dox_celltype="splicelib_hek_no_dox_rbp${rbp}"
    dox_celltype="splicelib_hek_dox_rbp${rbp}"

        echo "Running for RBP${rbp}: $dox_celltype and $no_dox_celltype"
        # Write the qsub command to the output file
        echo "qsub -v celltype1=$dox_celltype,celltype2=$no_dox_celltype /broad/dawnccle/melange/process_fastq/novaseq241106_pipeline/run_pairadise_pair_celltype_PSI.sh" >>"$output_file"
done



####### Process all cell line pairs #######
# Read the cell types into an array
readarray -t celltypes </broad/dawnccle/melange/process_fastq/novaseq241106_pipeline/novaseq241106_celllines.txt

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
        qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq/novaseq241106_pipeline/run_pairadise_pair_celltype_PSI.sh
    done
done

##### Process for 3ss also #####
# Read the cell types into an array
readarray -t celltypes </broad/dawnccle/melange/process_fastq/novaseq241106_pipeline/novaseq241106_celllines.txt

# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=()

# Specify the output file for the qsub commands
output_file="qsub_commands_3ss.txt"

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
        qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq/novaseq241106_pipeline/run_pairadise_pair_celltype_3ss.sh
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


###### Process for RBP type
readarray -t celltypes </broad/dawnccle/melange/data/RBP_shortlist_genes.txt

# Loop over each pair of cell types
for ((i = 0; i < ${#celltypes[@]}; i++)); do
        gene=${celltypes[i]}

        echo "Running for gene: $gene"
        # Write the qsub command to the output file
        qsub -v gene=$gene /broad/dawnccle/melange/process_fastq/missplicing/run_pairadise_pair_RBP_type_PSI.sh
        qsub -v gene=$gene /broad/dawnccle/melange/process_fastq/missplicing/run_pairadise_pair_RBP_type.sh
done

##### FOR RBP PSI #####
output_file="all_samples_RBP_PSI_combined_output.tsv"
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

##### FOR RBP 3ss #####
output_file="all_samples_RBP_3ss_combined_output.tsv"
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


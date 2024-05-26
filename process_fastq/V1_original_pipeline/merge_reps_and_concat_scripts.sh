#!/bin/bash

############################################################################################################
# Merge the samples based on the reps. 
############################################################################################################
# Extract unique sample-replicate combinations
combinations=$(find . -name "*.fastq.gz" | sed -E 's/.*\/([^-]+-rep[0-9]+)-.*/\1/' | sort -u)

for combo in $combinations; do
    echo $combo
    # Initialize arrays for R1 and R2 files for each combination
    R1_files=()
    R2_files=()

    # Find and add files to respective arrays
    for file in ${combo}-*; do
        if [[ $file =~ _R1_ ]]; then
            R1_files+=("$file")
        elif [[ $file =~ _R2_ ]]; then
            R2_files+=("$file")
        fi
    done

    # Merge R1 files for the current combination
    if [ ${#R1_files[@]} -gt 0 ]; then
        cat "${R1_files[@]}" > "${combo}_merged_R1_001.fastq.gz"
    fi

    # Merge R2 files for the current combination
    if [ ${#R2_files[@]} -gt 0 ]; then
        cat "${R2_files[@]}" > "${combo}_merged_R2_001.fastq.gz"
    fi
done



############################################################################################################
# Concatenate the dedup files.
############################################################################################################
# Output file
OUTPUT_FILE="concatenated_dedup_files.csv"

# Check if output file exists and remove it to start fresh
if [ -f "$OUTPUT_FILE" ]; then
    rm -f "$OUTPUT_FILE"
fi

# Variable to track if any files have been processed
FILES_PROCESSED=0

# Flag to indicate the first file
IS_FIRST_FILE=1

# Iterate over files matching the specified pattern
for file in *umi_dedup.csv; do
    echo $file
    # File exists and is a regular file, proceed
    FILES_PROCESSED=$((FILES_PROCESSED + 1))
    filename=$(basename "$file")

    # Remove the '_umi_dedup.csv' part from the filename
    filename_without_extension="${filename%_umi_dedup.csv}"

    # Use awk to append the modified filename as a new column to each line of the file
    # For the first file, include the header; for others, skip the first line
    if [ $IS_FIRST_FILE -eq 1 ]; then
        awk -v fname="$filename_without_extension" 'BEGIN{FS=OFS=","} {print $0, fname}' "$file" >>"$OUTPUT_FILE"
        IS_FIRST_FILE=0
    else
        awk -v fname="$filename_without_extension" 'BEGIN{FS=OFS=","} NR > 1 {print $0, fname}' "$file" >>"$OUTPUT_FILE"
    fi
done

if [ $FILES_PROCESSED -eq 0 ]; then
    echo "No files matched the pattern $PATTERN in directory $TARGET_DIR."
else
    echo "Concatenation complete. $FILES_PROCESSED files processed. Output saved to $OUTPUT_FILE."
fi

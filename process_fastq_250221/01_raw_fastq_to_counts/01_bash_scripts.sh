#!/bin/bash

# Merge the fastq files from 2 runs for the 230516 run. 
dir1="/broad/dawnccle/sequencing/merge1"
dir2="/broad/dawnccle/sequencing/merge2"
outdir="/broad/dawnccle/sequencing/merged_all"
mkdir -p "$outdir"

for read in R1 R2; do
  echo "Processing $read files..."

  for f2 in "$dir2"/*_${read}_001.fastq.gz; do
    [[ -e "$f2" ]] || continue

    base2=$(basename "$f2")
    prefix=$(echo "$base2" | sed -E "s/_merged_${read}_001\.fastq\.gz//")

    f1=$(find "$dir1" -type f -name "${prefix}-*_${read}_001.fastq.gz" | head -n1)

    if [[ -f "$f1" ]]; then
      outfile="$outdir/${prefix}_merged2_${read}_001.fastq.gz"
      echo "Concatenating: $(basename "$f2") + $(basename "$f1") -> $(basename "$outfile")"
      cat "$f2" "$f1" > "$outfile"
    else
      echo "WARNING: No match found in dir1 for $base2 using prefix: $prefix"
    fi
  done
done



dir1="/broad/dawnccle/sequencing/241106_Novaseq"
dir2="/broad/dawnccle/sequencing/250313_Novaseq"
outdir="/broad/dawnccle/sequencing/merged_all2"
mkdir -p "$outdir"

for read in R1 R2; do
  echo "Processing $read files..."

  for f2 in "$dir2"/splicelib_ZRSR2-rep*_${read}_001.fastq.gz; do
    [[ -e "$f2" ]] || continue

    base2=$(basename "$f2")
    # Extract sample ID: splicelib_ZRSR2-repX from splicelib_ZRSR2-rep1_S56_R1_001.fastq.gz
    sample=$(echo "$base2" | sed -E "s/_S[0-9]+_${read}_001\.fastq\.gz//")

    # Find matching file in dir1 with same sample ID and same read (ignore S numbers)
    f1=$(find "$dir1" -type f -name "${sample}_S*_${read}_001.fastq.gz" | head -n1)

    if [[ -f "$f1" ]]; then
      outfile="$outdir/${sample}_merged_${read}_001.fastq.gz"
      echo "Concatenating: $(basename "$f2") + $(basename "$f1") -> $(basename "$outfile")"
      cat "$f2" "$f1" > "$outfile"
    else
      echo "WARNING: No match found in dir1 for $base2 using sample: $sample"
    fi
  done
done

dir1="/broad/dawnccle/sequencing/241106_Novaseq"
dir2="/broad/dawnccle/sequencing/250313_Novaseq"
outdir="/broad/dawnccle/sequencing/merged_all2"
mkdir -p "$outdir"

for read in R1 R2; do
  echo "Processing $read files..."

  for f2 in "$dir2"/splicelib_U2AF1_WT-rep*_${read}_001.fastq.gz; do
    [[ -e "$f2" ]] || continue

    base2=$(basename "$f2")
    # Extract sample ID: splicelib_U2AF1_WT-repX from splicelib_U2AF1_WT-rep1_S56_R1_001.fastq.gz
    sample=$(echo "$base2" | sed -E "s/_S[0-9]+_${read}_001\.fastq\.gz//")

    # Find matching file in dir1 with same sample ID and same read (ignore S numbers)
    f1=$(find "$dir1" -type f -name "${sample}_S*_${read}_001.fastq.gz" | head -n1)

    if [[ -f "$f1" ]]; then
      outfile="$outdir/${sample}_merged_${read}_001.fastq.gz"
      echo "Concatenating: $(basename "$f2") + $(basename "$f1") -> $(basename "$outfile")"
      cat "$f2" "$f1" > "$outfile"
    else
      echo "WARNING: No match found in dir1 for $base2 using sample: $sample"
    fi
  done
done

dir1="/broad/dawnccle/sequencing/241106_Novaseq"
dir2="/broad/dawnccle/sequencing/250313_Novaseq"
outdir="/broad/dawnccle/sequencing/merged_all2"
mkdir -p "$outdir"

for read in R1 R2; do
  echo "Processing $read files..."

  for f2 in "$dir2"/splicelib_U2AF1_S34F-rep*_${read}_001.fastq.gz; do
    [[ -e "$f2" ]] || continue

    base2=$(basename "$f2")
    # Extract sample ID: splicelib_U2AF1_S34F-repX from splicelib_U2AF1_S34F-rep1_S56_R1_001.fastq.gz
    sample=$(echo "$base2" | sed -E "s/_S[0-9]+_${read}_001\.fastq\.gz//")

    # Find matching file in dir1 with same sample ID and same read (ignore S numbers)
    f1=$(find "$dir1" -type f -name "${sample}_S*_${read}_001.fastq.gz" | head -n1)

    if [[ -f "$f1" ]]; then
      outfile="$outdir/${sample}_merged_${read}_001.fastq.gz"
      echo "Concatenating: $(basename "$f2") + $(basename "$f1") -> $(basename "$outfile")"
      cat "$f2" "$f1" > "$outfile"
    else
      echo "WARNING: No match found in dir1 for $base2 using sample: $sample"
    fi
  done
done

dir1="/broad/dawnccle/sequencing/241106_Novaseq"
dir2="/broad/dawnccle/sequencing/250313_Novaseq"
outdir="/broad/dawnccle/sequencing/merged_all2"
mkdir -p "$outdir"

for read in R1 R2; do
  echo "Processing $read files..."

  for f2 in "$dir2"/splicelib_MEL202-rep*_${read}_001.fastq.gz; do
    [[ -e "$f2" ]] || continue

    base2=$(basename "$f2")
    # Extract sample ID: splicelib_MEL202-repX from splicelib_MEL202-rep1_S56_R1_001.fastq.gz
    sample=$(echo "$base2" | sed -E "s/_S[0-9]+_${read}_001\.fastq\.gz//")

    # Find matching file in dir1 with same sample ID and same read (ignore S numbers)
    f1=$(find "$dir1" -type f -name "${sample}_S*_${read}_001.fastq.gz" | head -n1)

    if [[ -f "$f1" ]]; then
      outfile="$outdir/${sample}_merged_${read}_001.fastq.gz"
      echo "Concatenating: $(basename "$f2") + $(basename "$f1") -> $(basename "$outfile")"
      cat "$f2" "$f1" > "$outfile"
    else
      echo "WARNING: No match found in dir1 for $base2 using sample: $sample"
    fi
  done
done

dir1="/broad/dawnccle/sequencing/241106_Novaseq"
dir2="/broad/dawnccle/sequencing/250313_Novaseq"
outdir="/broad/dawnccle/sequencing/merged_all2"
mkdir -p "$outdir"

for read in R1 R2; do
  echo "Processing $read files..."

  for f2 in "$dir2"/splicelib_sgCh3-rep*_${read}_001.fastq.gz; do
    [[ -e "$f2" ]] || continue

    base2=$(basename "$f2")
    # Extract sample ID: splicelib_sgCh3-repX from splicelib_sgCh3-rep1_S56_R1_001.fastq.gz
    sample=$(echo "$base2" | sed -E "s/_S[0-9]+_${read}_001\.fastq\.gz//")

    # Find matching file in dir1 with same sample ID and same read (ignore S numbers)
    f1=$(find "$dir1" -type f -name "${sample}_S*_${read}_001.fastq.gz" | head -n1)

    if [[ -f "$f1" ]]; then
      outfile="$outdir/${sample}_merged_${read}_001.fastq.gz"
      echo "Concatenating: $(basename "$f2") + $(basename "$f1") -> $(basename "$outfile")"
      cat "$f2" "$f1" > "$outfile"
    else
      echo "WARNING: No match found in dir1 for $base2 using sample: $sample"
    fi
  done
done
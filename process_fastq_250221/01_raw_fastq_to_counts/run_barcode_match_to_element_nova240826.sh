#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.
# This is good for jobs you need to run multiple times so you don't forget what it needs.

# Memory request for 4G
#$ -l h_vmem=30G

# Cores
#$ -pe smp 1
#$ -binding linear:1

# I like single output files
#$ -j y

# Not sure what this flag does
#$ -R y

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=12:00:00

# I don't like the top level of my homedir filling up.
#$ -o /broad/dawnccle/outputs/

# Job name
#$ -N Nova240826MatchBarcodeToElement

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse
source $HOME/.bioinfo
source activate python3.8

##################
### Run script ###
##################


# Run like:
# cd /broad/dawnccle/sequencing/240826_Novaseq
# for filename in *_R1_bc_extracted.fastq.gz; do
#     echo $filename
#     qsub -v FILENAME=$filename /broad/dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/run_barcode_match_to_element_nova240826.sh
# done


outdir=/broad/dawnccle/processed_data/reprocess_250221/nova240826
indir=/broad/dawnccle/sequencing/240826_Novaseq
mkdir -p $outdir

cd $indir

echo $FILENAME
BASENAME=${FILENAME%_merged_R1_001.fastq.gz}
FQ1=${BASENAME}_merged_R1_001.fastq.gz
FQ2=${BASENAME}_merged_R2_001.fastq.gz

BASENAME=${FILENAME%_R1_bc_extracted.fastq.gz}
FQ1=${BASENAME}_R1_bc_extracted.fastq.gz
FQ2=${BASENAME}_R2_bc_extracted.fastq.gz

python /broad/dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/MatchBarcodeToElementRNA_20250220.py \
    -1 ${indir}/${BASENAME}_R1_bc_extracted.fastq.gz \
    -l /broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv \
    -o $outdir

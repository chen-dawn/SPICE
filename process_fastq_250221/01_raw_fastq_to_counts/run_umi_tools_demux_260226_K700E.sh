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


# Not sure what this flag does
#$ -R y

# I like single output files
#$ -j y

# Job name
#$ -N K700E_UMI_TOOLS

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=24:00:00

# I don't like the top level of my homedir filling up.
#$ -o /broad/dawnccle/outputs/

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse
source $HOME/.bioinfo
source activate python3.8
# Use your dotkit
# use Python-3.4

##################
### Run script ###
##################

# Run like
# cd /broad/dawnccle/sequencing/260226_K700E
# for filename in lib47k*_R1_001.fastq.gz; do
#     echo $filename
#     qsub -v FILENAME=$filename /broad/dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/run_umi_tools_demux_260226_K700E.sh
# done

outdir=/broad/dawnccle/processed_data/reprocess_250221/260226_K700E
indir=/broad/dawnccle/sequencing/260226_K700E
mkdir -p $outdir

cd $indir

echo $FILENAME
BASENAME=${FILENAME%_R1_001.fastq.gz}
FQ1=${BASENAME}_R1_001.fastq.gz
FQ2=${BASENAME}_R2_001.fastq.gz


umi_tools extract --stdin $FQ1 --read2-in=${FQ2}  \
    --bc-pattern2="(?P<umi_1>.{10})(?P<discard_2>TTGCTAGGACCGGCCTTAAAGC){s<=2}(?P<cell_1>.{14})(?P<discard_3>CTACTTA){s<=1}.*" \
    --extract-method=regex \
    --whitelist /broad/thechenlab/Dawn/splicing/full_library_unique_exons_only/twist_47k_barcode_whitelist_reverse_complement.txt  \
    -L ${BASENAME}.log --stdout ${BASENAME}_R1_bc_extracted.fastq.gz --read2-out ${BASENAME}_R2_bc_extracted.fastq.gz \
    --filtered-out ${BASENAME}_R1_bc_extracted_failed.fastq.gz --filtered-out2 ${BASENAME}_R2_bc_extracted_failed.fastq.gz

python /broad/dawnccle/melange/process_fastq_250221/01_raw_fastq_to_counts/MatchBarcodeToElementRNA_20250220_original.py \
    -1 ${indir}/${BASENAME}_R1_bc_extracted.fastq.gz \
    -l /broad/dawnccle/melange/data/guide_library_cleaned/20240605_twist_library_v3_ID_barcode_ROUT_filtered.csv \
    -o $outdir
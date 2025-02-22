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
#$ -N UMI_TOOLS_AND_CALLING

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=24:00:00

# I don't like the top level of my homedir filling up.
#$ -o $HOME/outputs/

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

# Run like:
# while read p; do
#     echo $p
#     qsub -v FILENAME=$p /broad/dawnccle/melange/process_fastq/working_directory/umi_tools_then_splice_calling.sh
# done < merged_R1_filelist.txt

# For all R1 files in a folder:
# for filename in splicelib*_R1_001.fastq.gz; do
#     echo $filename
#     qsub -v FILENAME=$filename /broad/dawnccle/melange/process_fastq/working_directory/umi_tools_then_splice_calling.sh
# done

DIR=/broad/dawnccle/sequencing/240826_Novaseq
DIR=/broad/thechenlab/WalkUpSequencing/240826_SL-EXF_0201_A22KFHJLT3/Data/Intensities/BaseCalls/splicing
DIR=/broad/dawnccle/sequencing/241106_Novaseq/
cd $DIR

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


python /broad/dawnccle/melange/process_fastq/missplicing/MatchBarcodeToElementRNA_FINAL_VERSION_V5_240715.py \
    -1 ${DIR}/${BASENAME}_R1_bc_extracted.fastq.gz \
    -l /broad/dawnccle/melange/data/guide_library/20230130_twist_library_v3_ID_barcode_ROUT.csv \
    -o /broad/dawnccle/processed_data/novaseq_240826
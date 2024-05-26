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
#$ -N T-ZUMI_TOOLS

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=72:00:00

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

cd /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/merged_fastqs/

for FILENAME in [T-Z]*_merged_R1_001.fastq.gz; do
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

    # FQ1=${BASENAME}_R1_bc_extracted.fastq.gz 
    # FQ2=${BASENAME}_R2_bc_extracted.fastq.gz 
    # python /broad/thechenlab/Dawn/splicing/library_sequencing_47k/MatchBarcodeToElementRNA_umi_tools_extracted_47k.py \
    #     -1 /broad/thechenlab/NextSeqOutput_Xerneas/230308_VL00297_152_AACG5CGM5/Data/Intensities/BaseCalls/${FQ1} \
    #     -2 /broad/thechenlab/NextSeqOutput_Xerneas/230308_VL00297_152_AACG5CGM5/Data/Intensities/BaseCalls/${FQ2} \
    #     -l /broad/thechenlab/Dawn/splicing/library_sequencing_47k/20230130_twist_library_v3_ID_barcode_ROUT.csv \
    #     -r /broad/thechenlab/Dawn/splicing/library_sequencing_47k/RNA_ref/WEAK_47k_reference.fasta \
    #     -o /broad/thechenlab/Dawn/splicing/library_sequencing_47k/230521_RNA_celltypes

done


# count number of reads in each fastq file.
# rm fastq_read_counts.txt
# for i in `ls *.fastq.gz`; do 
# echo $i
# echo -en $i >> fastq_read_counts.txt
# echo -en "\t" >> fastq_read_counts.txt
# echo $(zcat ${i} | wc -l)/4|bc >> fastq_read_counts.txt

# done

# Count the number of unique cell barcodes in each fastq file.
# for i in `ls *R2_bc_extracted.fastq.gz`; do 
# echo $i
# BASENAME=${i%_R2_bc_extracted.fastq.gz}
# zgrep ^@V $i | awk -F'[ _]' '{print $2"_"$3}' | sort | uniq -c > ${BASENAME}_cell_barcode_umi_counts.txt
# zgrep ^@V $i | awk -F'[ _]' '{print $2}' | sort | uniq -c > ${BASENAME}_cell_barcode_counts.txt
# done

# Parallel version of that
# count_barcode_umi_task(){
# FILENAME=$1
# echo $FILENAME
# BASENAME=${FILENAME%_R2_bc_extracted.fastq.gz}
# zgrep ^@V $1 | awk -F'[ _]' '{print $2"_"$3}' | sort | uniq -c > ${BASENAME}_cell_barcode_umi_counts.txt
# zgrep ^@V $1 | awk -F'[ _]' '{print $2}' | sort | uniq -c > ${BASENAME}_cell_barcode_counts.txt
# }

# N=8
# (
# for FILE in *R2_bc_extracted.fastq.gz; do
# ((i=i%N)); ((i++==0)) && wait
# count_barcode_umi_task "$FILE" &
# done
# )

# # Count the number of unique cell barcodes in each fastq file.
# for i in `ls *cell_barcode_counts.txt`; do 
# wc -l $i >> merged_cell_barcode_counts.txt
# done
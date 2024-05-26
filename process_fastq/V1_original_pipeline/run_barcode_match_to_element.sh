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
#$ -l h_rt=24:00:00

# I don't like the top level of my homedir filling up.
#$ -o $HOME/outputs/

# Job name
#$ -N MatchBarcodeToElement

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
# while read p; do
#     echo $p
#     qsub -v filename=$p run_barcode_match_to_element.sh
# done < all_fastq_file_names.txt

# cd /broad/thechenlab/Dawn/Helicase/221215ExomeSeq

filebase=${filename%.bam}
filebase=${filename%_R1_bc_extracted.fastq.gz}

python /broad/thechenlab/Dawn/splicing/library_sequencing_47k/MatchBarcodeToElementRNA_umi_tools_extracted_47k_R1_VCP_230524.py \
    -1 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/${filename} \
    -2 /broad/dawnccle/230516_SL-EXC_0008_B2235L7LT3/Data/Intensities/BaseCalls/${filebase}_R2_bc_extracted.fastq.gz \
    -l /broad/thechenlab/Dawn/splicing/library_sequencing_47k/20230130_twist_library_v3_ID_barcode_ROUT.csv \
    -r /broad/thechenlab/Dawn/splicing/library_sequencing_47k/RNA_ref/WEAK_47k_reference.fasta \
    -o /broad/dawnccle/processed_data/celltype47_230524
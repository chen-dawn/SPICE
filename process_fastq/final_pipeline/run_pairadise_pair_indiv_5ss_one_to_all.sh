#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.
# This is good for jobs you need to run multiple times so you don't forget what it needs.

# Memory request for 4G
#$ -l h_vmem=4G

# Cores
#$ -pe smp 8
#$ -binding linear:8

# I like single output files
#$ -j y

# Not sure what this flag does
#$ -R y

# Runtime request.  Usually 30 minutes is plenty for me and helps me get backfilled into reserved slots.
#$ -l h_rt=24:00:00

# I don't like the top level of my homedir filling up.
#$ -o $HOME/ccle/outputs

# Job name
#$ -N RMATS_5ss_FINAL

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse
source $HOME/.bioinfo
reuse R-4.1

##################
### Run script ###
##################

# celltype1=KMRC1
# celltype2=HCC38

echo $celltype1

output_dir=/broad/dawnccle/processed_data/latest/rmats_stat_out_indiv_5ss/$celltype1
output_file=/broad/dawnccle/processed_data/latest/rmats_stat_out_indiv_5ss/$celltype1/$celltype1\_formatted_df.tsv

mkdir -p $output_dir
Rscript /broad/dawnccle/melange/process_fastq/final_pipeline/pairadise_indiv_5ss_one_to_all.R $celltype1 $output_file

source activate rmats2.7

python /broad/dawnccle/rMATS-STAT/rMATS_unpaired.py $output_file $output_dir 8 0.1
python /broad/dawnccle/rMATS-STAT/FDR.py $output_dir/rMATS_Result_P.txt $output_dir/$celltype1\_rMATS_Result_P.FDR.txt
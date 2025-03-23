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
#$ -o /broad/dawnccle/outputs/

# Job name
#$ -N PSI_RBP_RMATS_STAT

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

# Run like this:
# qsub_commands_file=/broad/dawnccle/melange/process_fastq_250221/06_analyze_RBPs/qsub_commands_RBP_shortlist_genes.txt
# > $qsub_commands_file
# while read gene; do
#     echo $gene
#     echo "qsub -v gene=$gene /broad/dawnccle/melange/process_fastq_250221/06_analyze_RBPs/run_pairadise_pair_RBP_type_PSI.sh" >> $qsub_commands_file
# done < /broad/dawnccle/processed_data/reprocess_250221/rbp_gex_outputs/RBP_shortlist_genes.txt

echo $gene
output_dir=/broad/dawnccle/processed_data/reprocess_250221/pairadise_RBP_PSI/RBP_type_PSI_broad/$gene
output_file=/broad/dawnccle/processed_data/reprocess_250221/pairadise_RBP_PSI/RBP_type_PSI_broad/$gene/$gene\_formatted_df.tsv

mkdir -p $output_dir
Rscript /broad/dawnccle/melange/process_fastq_250221/06_analyze_RBPs/pairadise_PSI_RBP_type.R $gene $output_file
source activate rmats2.7

python /broad/dawnccle/rMATS-STAT/rMATS_unpaired.py $output_file $output_dir 8 0.1
python /broad/dawnccle/rMATS-STAT/FDR.py $output_dir/rMATS_Result_P.txt $output_dir/$gene\_rMATS_Result_P.FDR.txt
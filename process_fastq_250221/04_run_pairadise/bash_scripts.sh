
######################
# This is for OEx pairs
######################
# Compare the dox vs no dox only
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_PSI_OEx.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Loop over RBP identifiers (1 to 12)
for rbp in {1..12}; do
    # Define the patterns for dox and no_dox for the current RBP
    no_dox_celltype="splicelib_hek_no_dox_rbp${rbp}"
    dox_celltype="splicelib_hek_dox_rbp${rbp}"

        echo "Running for RBP${rbp}: $dox_celltype and $no_dox_celltype"
        # Write the qsub command to the output file
        echo "qsub -v celltype1=$dox_celltype,celltype2=$no_dox_celltype /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_OEx_PSI.sh" >>"$output_file"
done

######################
# This is for MUT pairs
######################
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_PSI_MUT.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Get all cell types in the MUT directory
celltypes=($(ls /broad/dawnccle/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/MUT/*.tsv | xargs -n 1 basename | sed -E 's/-rep[0-9]+.*//' | sort -u))
# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=()

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
        echo "qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_MUT_PSI.sh" >>"$output_file"
    done
done

######################
# This is for MUT2 pairs
######################
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_PSI_MUT2.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Get all cell types in the MUT directory
celltypes=($(ls /broad/dawnccle/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/MUT2/*.tsv | xargs -n 1 basename | sed -E 's/-rep[0-9]+.*//' | sort -u))
# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=()

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
        echo "qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_MUT2_PSI.sh" >>"$output_file"
    done
done


######################
# This is for WT pairs one to one
######################
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_PSI_WT_one_to_one.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Get all cell types in the WT directory
celltypes=($(ls /broad/dawnccle/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/WT/*.tsv | xargs -n 1 basename | sed -E 's/-rep[0-9]+.*//' | sort -u))
# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=() 

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
        echo "qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_WT_PSI.sh" >>"$output_file"
    done
done

######################
# This is for WT pairs one to all
######################
# Specify the output file for the qsub commands. Make full path
output_file="/broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/qsub_commands_PSI_WT_one_to_all.txt"
# Ensure the output file is empty at the start
>"$output_file"

# Get all cell types in the WT directory
celltypes=($(ls /broad/dawnccle/processed_data/reprocess_250221/count_normalized_chimeric_rate_considering_included/WT/*.tsv | xargs -n 1 basename | sed -E 's/-rep[0-9]+.*//' | sort -u))
# Create a set to track processed pairs
declare -A processed_pairs
processed_pairs=()

#  Loop over each cell type
for celltype in "${celltypes[@]}"; do
    echo "Running for celltype: $celltype"
    # Write the qsub command to the output file
    echo "qsub -v celltype1=$celltype,celltype2=all /broad/dawnccle/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_WT_PSI.sh" >>"$output_file"
done


######################
# This is for MUT2 pairs but on workstation
######################

SCRIPT="/mnt/dawnccle2/melange/process_fastq_250221/04_run_pairadise/run_pairadise_pair_MUT2_PSI_workstation.sh"

declare -a PAIRS=(
"CH3-1_A1 CH3-1_A2"
"CH3-1_A1 FUBP1_B12"
"CH3-1_A1 FUBP1_C5"
"CH3-1_A1 RBM10_C8"
"CH3-1_A1 RBM10_G4"
"CH3-1_A1 RBM5_A2"
"CH3-1_A1 RBM5_A3"
"CH3-1_A1 ZRSR2_F8"
"CH3-1_A1 ZRSR2_G9"
"CH3-1_A1 splicelib_sgCh3"
"CH3-1_A2 FUBP1_B12"
"CH3-1_A2 FUBP1_C5"
"CH3-1_A2 RBM10_C8"
"CH3-1_A2 RBM10_G4"
"CH3-1_A2 RBM5_A2"
"CH3-1_A2 RBM5_A3"
"CH3-1_A2 ZRSR2_F8"
"CH3-1_A2 ZRSR2_G9"
"CH3-1_A2 splicelib_sgCh3"
"FUBP1_B12 FUBP1_C5"
"FUBP1_B12 RBM10_C8"
"FUBP1_B12 RBM10_G4"
"FUBP1_B12 RBM5_A2"
"FUBP1_B12 RBM5_A3"
"FUBP1_B12 ZRSR2_F8"
"FUBP1_B12 ZRSR2_G9"
"FUBP1_C5 RBM10_C8"
"FUBP1_C5 RBM10_G4"
"FUBP1_C5 RBM5_A2"
"FUBP1_C5 RBM5_A3"
"FUBP1_C5 ZRSR2_F8"
"FUBP1_C5 ZRSR2_G9"
"K562WT K562K700E"
"RBM10_C8 RBM10_G4"
"RBM10_C8 RBM5_A2"
"RBM10_C8 RBM5_A3"
"RBM10_C8 ZRSR2_F8"
"RBM10_C8 ZRSR2_G9"
"RBM10_G4 RBM5_A2"
"RBM10_G4 RBM5_A3"
"RBM10_G4 ZRSR2_F8"
"RBM10_G4 ZRSR2_G9"
"RBM5_A2 RBM5_A3"
"RBM5_A2 ZRSR2_F8"
"RBM5_A2 ZRSR2_G9"
"RBM5_A3 ZRSR2_F8"
"RBM5_A3 ZRSR2_G9"
"splicelib_sgCh3 splicelib_ZRSR2"
"splicelib_U2AF1_WT splicelib_U2AF1_S34F"
"splicelib_ZRSR2 ZRSR2_F8"
"splicelib_ZRSR2 ZRSR2_G9"
"ZRSR2_F8 ZRSR2_G9"
)

for pair in "${PAIRS[@]}"; do
  celltype1=$(echo $pair | cut -d' ' -f1)
  celltype2=$(echo $pair | cut -d' ' -f2)
  echo "Running: $celltype1 vs $celltype2"
  celltype1="$celltype1" celltype2="$celltype2" $SCRIPT
done


# Run the RBP types on workstation.
SCRIPT="/mnt/dawnccle2/melange/process_fastq_250221/06_analyze_RBPs/run_pairadise_pair_RBP_type_PSI_workstation.sh"
declare -a GENES=(
  PIH1D2 ANK3 MAT1A MMAA PABPC3 RBM46 BOLL ZFP36L2 DDX4 JAKMIP1 SRSF12 SORBS2 SAMSN1 PPARGC1B AFF2 TDRD9 F11R NPM2 EPB41
  IGF2BP1 CELF3 GNE AZGP1 CELF5 DCD ZNF385A ELAVL4 HENMT1 TDRD5 S100A9 TDRD10 EIF4E3 IFI16 RBM47 ZC3H12A ADAD1 TLR3 TERT
  ZCCHC24 RPL10L NDRG2 RNASE7 FBN1 RBPMS2 RBFOX3 XIRP1 MTCL1 ROR2 PCSK9 NR0B1 GP2 RNASE2 RNASE3 RNASE6 RBMY1F SMAD1 RBMXL2
  PLA2G1B NCBP2L SNTB1 ISG20 CSDC2 ZMAT3 PURG RNASE8 RNASE11 APOBEC4 TDRD12 PABPC5 RBMXL3 EIF4E1B BASP1 RBM44 SRRM3 DDX10
  ZC3H12D ERN1 KCTD12 RPP25 HNRNPCL1 NSUN7 MKRN3 PLD6 APOBEC3B NLRP11 TDRD6 CFAP65 DDX60L IBA57 RNASE10 CADM1 MEX3B PCBP3
  PABPC1L2B PIWIL3 RALYL DDX53 RBM11 PDIA2 ADARB2 IFIT1 PABPC1L2A MAPT DAZ3 NANOS3 LIN28B COL14A1 CPSF4L ZCCHC13 DAZ1
  ZC3H6 NANOS2 FBLL1 NANOS1 ADAT2 GSPT2 TDRD7 S100A4 ELAVL3 MYO18A LAMA2 TLR7 ZNF239 PIWIL2 ADARB1 MYO5A HMGN5 EIF1AY
  EPS8L3 RBM20 IFIT1B HSPA1A POU5F1 RANBP17 PABPN1L DAZ4 DAZ2 NYNRIN RNASE13 TRIM71 ANG CPEB1 ARHGEF28 TDRD15 PAPOLB RBMY1J
  POLR2J2 PATL2 RBMY1A1 APOBEC3G L1TD1 PWP2 RBMY1E RBMY1B RBMY1D ZCCHC3 RBM14-RBM4 EIF5AL1 UTP14C PABPC4L MEX3A FDXACB1 DND1
  KHDC1L RNASE4 EPPK1 PCDHGA9 TEX13A CALR3 FBXO17 NXF2 NXF2B NBPF10 C2orf15 LENG9 RDM1 PCDH20 RPS4Y2
)

for gene in "${GENES[@]}"; do
  echo "Running: $gene"
  gene="$gene" $SCRIPT
done

# Run the RBP types manually on uger1.
SCRIPT="/broad/dawnccle/melange/process_fastq_250221/06_analyze_RBPs/run_pairadise_pair_RBP_type_PSI.sh"
declare -a GENES=(
  YBX2 DCN SAMD4A RNASET2 CELF2 DDX3Y HEATR6 SIDT1 MOV10L1 MBNL3 MAP2 RBFOX1 TNRC6C TNS1 DDX1 RIMS1 DDX43 AFP TDRD3 APOB DNMT3B
  P2RX7 OAS1 DAZL TDRD1 SORBS1 DSP APOBEC3H CHGA EEF1A2 SAMHD1 CELF4 NXT2 TLR8 LUZP4 ZC3H12B CORO1A NDRG4 ESRP2 SCG3 ESRP1
  NDRG1 NOVA2 ZFR2 MYH14 OGN DNM1 ELAVL2 DDX58 CPEB3 HOXB6 DHX58 PTGES3L-AARSD1 ALDOC PPARGC1A DDX25 CRYAB OAS3 OAS2 ENDOU
  )

for gene in "${GENES[@]}"; do
  echo "Running: $gene"
  gene="$gene" $SCRIPT
done

# Run the RBP types manually on uger2.
SCRIPT="/broad/dawnccle/melange/process_fastq_250221/06_analyze_RBPs/run_pairadise_pair_RBP_type_PSI.sh"
declare -a GENES=(
  APOBEC1 RBM24 KHDRBS2 IQCG IFIH1 FN1 SWT1 LEPR WARS2 ALDH6A1 IFIT3 IFIT2 ENOX1 SMAD9 UTP20 BICC1 LRP1 C4BPA PAIP2B APOBEC2
  ATXN1 PIWIL1 AGO3 NXF5 VIL1 ZFP36 APOBEC3F RNASE1 RPS4Y1 CNN1 BST2 HELZ2 ASS1 DUSP9 SYNE1 MAP1B KHDRBS3 LIN28A LGALS3
  DMGDH RNF17 CTIF ERN2 PIWIL4 HOOK1 FADS2 DZIP1 MSI1 OASL TES KHDC1 RNASEL CPEB2 DDX60 SMAD6 SPTBN5 HERC5 SRRM4 RNF113B NOVA1 CELF6
  ADAD2 RPL3L SLC16A3 MAEL CGN DQX1 AFF3 RBMS3 SLC25A48 CPLX2 PNLDC1 ZC3HAV1L NXF3 CDKN2A A1CF GAS2 ZC3H12C TAGLN PDCD4
  )

for gene in "${GENES[@]}"; do
  echo "Running: $gene"
  gene="$gene" $SCRIPT
done


# ##### Process for 3ss also #####
# # Read the cell types into an array
# readarray -t celltypes </broad/dawnccle/melange/process_fastq/novaseq241106_pipeline/novaseq241106_celllines.txt

# # Create a set to track processed pairs
# declare -A processed_pairs
# processed_pairs=()

# # Specify the output file for the qsub commands
# output_file="qsub_commands_3ss.txt"

# # Ensure the output file is empty at the start
# >"$output_file"

# # Loop over each pair of cell types
# for ((i = 0; i < ${#celltypes[@]}; i++)); do
#     for ((j = i + 1; j < ${#celltypes[@]}; j++)); do
#         celltype1=${celltypes[i]}
#         celltype2=${celltypes[j]}

#         # Skip if celltype1 is the same as celltype2
#         if [[ "$celltype1" == "$celltype2" ]]; then
#             continue
#         fi

#         # Create a sorted key to ensure pairs are unique
#         pair_key=$(echo "$celltype1 $celltype2" | tr " " "\n" | sort | tr "\n" "_")

#         # Skip if the pair has already been processed
#         if [[ -n "${processed_pairs[$pair_key]}" ]]; then
#             continue
#         fi

#         # Mark the pair as processed
#         processed_pairs[$pair_key]=1

#         echo "Running for pair: $celltype1 and $celltype2"
#         # Write the qsub command to the output file
#         qsub -v celltype1=$celltype1,celltype2=$celltype2 /broad/dawnccle/melange/process_fastq/novaseq241106_pipeline/run_pairadise_pair_celltype_3ss.sh
#     done
# done



# Merge the sequences for FDR for MUT
output_file="MUT_PSI_combined_output_indiv.tsv"
first_file=true

for folder in MUT_*/; do
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

# Merge the sequences for FDR for OEx
output_file="OEx_PSI_combined_output_indiv.tsv"
first_file=true

for folder in OEx_*/; do
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


# Merge the sequences for FDR for WT.
output_file="WT_PSI_combined_output_indiv.tsv"
first_file=true

for folder in WT_*/; do
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

# Merge the sequences for FDR for MUT2
output_file="MUT2_PSI_combined_output_indiv.tsv"
first_file=true

for folder in MUT2_*/; do
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
    

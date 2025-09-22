#!/bin/bash
for FILENAME in satmutv2*_R1_bc_extracted.fastq.gz; do
    echo $FILENAME
    BASENAME=${FILENAME%_R1_bc_extracted.fastq.gz}
    FQ1=${BASENAME}_R1_bc_extracted.fastq.gz
    FQ2=${BASENAME}_R2_bc_extracted.fastq.gz

    python /broad/thechenlab/Maile/splicing/sat_mut_scripts/MatchBarcodeToElementRNA_250703_no_mappy.py \
        -1 /broad/dawnccle/sequencing/250701_NextSeq_satmutv2/${FQ1} \
        -l /broad/thechenlab/Maile/splicing/mutagenesis_library/20250604_twist_library_mutagenesis_barcodes.csv \
        -o /broad/dawnccle/processed_data/satmutv2_dawn
done

cd /broad/dawnccle/sequencing/241106_Novaseq/satmut_unfiltered/
for FILENAME in satmut*_R1_bc_extracted.fastq.gz; do
    echo $FILENAME
    BASENAME=${FILENAME%_R1_bc_extracted.fastq.gz}
    FQ1=${BASENAME}_R1_bc_extracted.fastq.gz
    FQ2=${BASENAME}_R2_bc_extracted.fastq.gz

    python /broad/thechenlab/Maile/splicing/sat_mut_scripts/MatchBarcodeToElementRNA_250703_no_mappy.py \
        -1 /broad/dawnccle/sequencing/241106_Novaseq/satmut_unfiltered/${FQ1} \
        -l /broad/thechenlab/Maile/splicing/mutagenesis_library/20240901_twist_library_mutagenesis_barcodes.csv \
        -o /broad/dawnccle/processed_data/satmutv1_dawn
done
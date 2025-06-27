#! /bin/bash

nextflow run nf-core/rnasplice \
  --input /mnt/dawnccle2/melange/process_fastq_250221/09_analyze_RBP_overexpression/samplesheet.csv \
  --contrasts /mnt/dawnccle2/melange/process_fastq_250221/09_analyze_RBP_overexpression/contrastsheet.csv \
  --fasta /mnt/dawnccle2/reference/gencode/GRCh38.primary_assembly.genome.fa \
  --gtf /mnt/dawnccle2/reference/gencode/gencode.v48.primary_assembly.annotation.gtf \
  --aligner star \
  --outdir /mnt/dawnccle2/nf_rnasplice_OEx \
  -r 1.0.4 \
  -profile docker \
   -work-dir /data/dc_tmp \
  --max_cpus 64 \
  --max_memory '256.GB' \
  -process.executor local 

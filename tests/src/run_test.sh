#!/bin/bash

proj_dir="${PWD%/*}"

nextflow ../../main.nf \
    --publish_dir ${proj_dir}/results/ill_test \
    --genome ${proj_dir}/data/fortune_primary_v1.4.0.fasta.masked \
    --fastq_pe /pickett_flora/projects/citrus/raw_data/RNASEQ_c_reticulata_PRJNA999253 \
    --skip_qc true \
    --skip_trim true \
    -profile local,four

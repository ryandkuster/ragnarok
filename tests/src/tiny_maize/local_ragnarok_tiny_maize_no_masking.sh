#!/bin/bash

source config
echo -e "${ANALYSIS}\n"
echo "project  : ${PROJ}"
echo "data     : ${DATA}"
echo "src      : ${SRC}"
echo "results  : ${RESULTS}"
echo "derived  : ${DERIVED}"

cds=$DATA/tiny_maize/subset_Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.fna
genome=$DATA/tiny_maize/chr1_171332842_173553286.fa
protein=$DATA/tiny_maize/subset_Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.faa
ill=$DATA/tiny_maize/fastq/

RAG="${PWD%/*/*/*}"

nextflow  ../../../main.nf \
    --design          $RAG/assets/mikado_conf.tsv \
    --scoring         $RAG/assets/plant_microexon_filter.yaml \
    --publish_dir     $RESULTS/local_run/no_mask \
    --genome          $genome \
    --protein         $protein \
    --ill             $ill \
    --scoring         $DATA/db/plant.yaml \
    --homology        $DATA/db/uniprotkb_taxonomy_id_33090_AND_reviewe_2025_03_07.fasta \
    --skip_qc         false \
    --skip_trim       false \
    --nlrs            false \
    --perform_masking false \
    --cds             $cds \
    --entap_db        $DATA/db/entap_db \
    --final_prefix    no_mask \
    -profile          local,eight \
    -resume



#!/bin/bash
#SBATCH -J nf_annot
#SBATCH -A acf-utk0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=0-12:00:00
#SBATCH --error=logs/job.e%J
#SBATCH --output=logs/job.o%J

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

nextflow  ../../../main.nf \
    --design          ~/nextflow/ragnarok/assets/mikado_conf.tsv \
    --publish_dir     $RESULTS/slurm_run \
    --genome          $genome \
    --protein         $protein \
    --ill             $ill \
    --scoring         $DATA/db/plant.yaml \
    --homology        $DATA/db/uniprotkb_taxonomy_id_33090_AND_reviewe_2025_03_07.fasta \
    --skip_qc         false \
    --skip_trim       false \
    --nlrs            false \
    --perform_masking true \
    --cds             $cds \
    -profile          slurm,custom \
    -resume



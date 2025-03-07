# Ragnarok

```
################################################################################
#                                   RAGNAROK!                                  #
#                         RApid Genome anNotAtion ROcKs!                       #
#                          Chris Gottschalk 12/6/2024                          #
#                                    v1.1                                      #
################################################################################
```

## quick start guide

### requirements

nextflow (22.10.4+)  
apptainer (1.1.8+)  
path to helixer model https://zenodo.org/records/10836346

### getting started

### running on a local server

First, do you have a working nextflow/apptainer version?

```
nextflow -version
apptainer --version
```

```
nextflow <path to your pipeline>/ragnarok/main.nf \
    --publish_dir <your desired results directory> \
    --genome <path to reference genome> \
    --skip_qc false \
    --skip_trim false \
    --skip_mark_dupe true \
    -profile local,four \
    -resume
```

Note the profile here is set up for use on a local server, but will likely require modification for your job. The `four` local profile is set up to use approximately 4 cpus maximum. Other presets exist in `conf/local.conf` and you can create your own by copying those examples.


### running on a slurm server

First, do you have a working nextflow/apptainer version?

```
nextflow -version
apptainer --version
```

Below is a sample sbatch script to run the pipeline. You'll need to replace the values with those that make sense for your use case.

```
#!/bin/bash
#SBATCH -J ragnarok
#SBATCH -A acf-utk0032
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=5-00:00:00
#SBATCH --error=job.e%J
#SBATCH --output=job.o%J

module load nextflow/23.10.0

export NXF_OPTS="-Xms500M -Xmx2G"
export NXF_ANSI_LOG=false

nextflow  <path to your pipeline>/vary_cool/main.nf \
    --publish_dir <your desired results directory> \
    --genome <path to reference genome> \
    --skip_qc false \
    --skip_trim false \
    -profile slurm,custom \
    -resume
```

Based on your qos/partition, you may want to modify the `conf/slurm.config` and `conf/slurm_custom.config` files to handle your dataset.

#### a few notes on slurm qos and partitions

This pipeline currently relies on four qos/partition configurations on the UTK ISAAC-NG system.

- short : maximum 3 hours, 12 jobs submitted
- campus : maximum 1 day, 94 jobs submitted
- long : maximum 6 days, 14 jobs submitted
- gpu : maximum 6 days, 14 jobs submitted

For each of the three labels found in `slurm.config` maxForks can be adjusted to qos that work on other sytems, and the imported clusterOptions for each can be updated with specific qos/partition/account information found in `slurm_custom.config`.

To check the limits for a given qos, replace short with any qos:

```
sacctmgr show qos where name=short
```

...or

```
scontrol show partition short
```

## development planning:

### desired input considerations

### necessary steps  
- [x] alignment 
  - [x] alignment - STAR
  - [ ] alignment - long reads
- [x] helixer
- [x] stringtie
- [x] gffread
- [x] miniprot
- [x] mikado2
- [x] busco
- [x] handle additional gff input paths
- [ ] determine steps where copying to publish_dir is needed

### optional steps  
- [x] QC
- [x] trimming
- [ ] FindPlantNLRs annotation
- [ ] EDTA masking
- [x] allow non land_plant model for Helixer (opened issue)

### obstacles
- [ ] mikado2 quay container is broken
- [ ] edta run fails on citrus genome

```mermaid
flowchart TB
    subgraph " "
    subgraph params
    v25["lineage"]
    v15["read_type"]
    v38["scoring"]
    v39["homology"]
    v11["masked"]
    v2["skip_qc"]
    v6["minimum_length"]
    v30["skip_hx"]
    v12["genome"]
    v13["cds"]
    v0["fastq_pe"]
    v23["protein"]
    v31["subseq_len"]
    v34["design"]
    v5["skip_trim"]
    v19["skip_st"]
    end
    v3([FASTQC_RAW])
    v4([MULTIQC_RAW])
    v7([FASTP_ADAPTERS])
    v9([FASTQC_TRIM])
    v10([MULTIQC_TRIM])
    v14([EDTA])
    v16([STAR_INDEX_NA])
    v17([STAR_MAP])
    v18([SAM_SORT])
    v20([STRINGTIE])
    v21([GFFREAD])
    v22([TRANSDECODER])
    v24([MINIPROT])
    v26([HELIXER_DB])
    v32([HELIXER])
    v35([PARSE_INPUT])
    v40([MIKADO_CONF])
    v41([TRANSDECODER_ORF])
    v42([DIAMOND])
    v43([THE_GRANDMASTER])
    v44([GFFREAD_FINAL])
    v45([BUSCO])
    v46([COMPLEASM])
    v0 --> v3
    v3 --> v4
    v0 --> v7
    v6 --> v7
    v7 --> v9
    v9 --> v10
    v12 --> v14
    v13 --> v14
    v12 --> v16
    v16 --> v17
    v7 --> v17
    v17 --> v18
    v18 --> v20
    v20 --> v21
    v12 --> v21
    v20 --> v22
    v12 --> v22
    v22 --> v24
    v23 --> v24
    v12 --> v24
    v25 --> v26
    v26 --> v32
    v12 --> v32
    v31 --> v32
    v34 --> v35
    v19 --> v35
    v30 --> v35
    v32 --> v40
    v35 --> v40
    v21 --> v40
    v38 --> v40
    v22 --> v40
    v39 --> v40
    v24 --> v40
    v12 --> v40
    v40 --> v41
    v39 --> v42
    v40 --> v42
    v39 --> v43
    v40 --> v43
    v41 --> v43
    v42 --> v43
    v12 --> v43
    v43 --> v44
    v12 --> v44
    v44 --> v45
    v44 --> v46
    end
```
# Ragnarok

```
################################################################################
#                                   RAGNAROK!                                  #
#                         RApid Genome anNotAtion ROcKs!                       #
#                          Chris Gottschalk 12/6/2024                          #
#                                    v1.1                                      #
################################################################################
```

# requirements

## software
nextflow (22.10.4+)  
apptainer (1.1.8+)

## required files
**genome** : Provide your assembly (ideally with *simple* names if using EDTA masking option).  
**short/long reads** : any combination of short/long reads can be used for input.  
**helixer** : Helixer model https://zenodo.org/records/10836346 (land_plants default).  
**miniprot** : Protein file for alignment (e.g., closest ref species).  
**mikado2** : Mikado protein homology file (e.g., [uniprot 33090 for viridiplantae](https://www.uniprot.org/uniprotkb?query=viridiplantae&facets=reviewed%3Atrue)).  
**mikado2** : Mikado scoring file (e.g., [plant.yaml](https://github.com/EI-CoreBioinformatics/mikado/tree/master/Mikado/configuration/scoring_files)).  
**mikado2** : Mikado configuration table file (see Mikado [documentation](https://mikado.readthedocs.io/en/stable/Tutorial/)).  

The mikado2 stage of the pipeline requires a configuration table to weigh the input gene models and give model priority. For this `--design` input, Ragnarok has the mandatory fields (`hx`, `st`, `mp`, `tr`) for the helixer, stringtie, miniprot, and transdecoder models produced along the way. Users should leave the file field blank if they are to be performed in the pipeline, but all other fields should be present. 

The fields in this file are `file location`, `alias`, `strand-specific`, `sample-score`, `reference`, and `exclude redundant models`. Mikado [documentation](https://mikado.readthedocs.io/en/stable/Tutorial/) provides further information on the choice of values in this table.

Example tsv configuration (assets/mikado_conf.tsv):

```
	hx	True		False	False
	st	True	1	False	True
	tr	False	-0.5	False	False
	mp	True	1	False	False
```

> [!NOTE]
> _The first field is intentionally missing as Ragnarok will produce these outputs. _

Example tsv configuration for `--nlrs true` (assets/mikado_nlr_conf.tsv):

```
	hx	True		False	False
	st	True	1	False	True
	tr	False	-0.5	False	False
	mp	True	1	False	False
	nlr	True	1	False	False
```

> [!NOTE]
> _The first field is intentionally missing as Ragnarok will produce these outputs. _

Ragnarok also allows for any number of *existing* input annotations (gff3) to be input as additional models into the mikado2 stage of processing.

Hypothetical tsv configuration including combination of Ragnarok (empty file fields) and existing models:

```
	hx	True		False	False
	st	True	1	False	True
	tr	False	-0.5	False	False
	mp	True	1	False	False
cufflinks.gtf	cuff	True		False	False
trinity.gff3	tr	False	-0.5	False	False
reference.gff3	at	True	5	True	False
```

> [!NOTE]
> _The filepath to existing gffs will need to be provided. _


## optional files

**edta** : cds file for your species (used for masking)  
**findplantnlrs** : interproscan ([64-bit download](https://www.ebi.ac.uk/interpro/download/InterProScan/))  
**findplantnlrs** : genemark configured for container (see assets/genemark_setup.sh)

# getting started

## running on a local server

First, do you have a working nextflow/apptainer version?

```
nextflow -version
apptainer --version
```

Below is a sample script to run the pipeline. You'll need to replace the `<>` values with those that make sense for your use case.

```
nextflow  ~/nextflow/ragnarok/main.nf \
    --publish_dir < path to results location > \
    --genome      < path to genome in fasta > \
    --cds         < path to cds fasta file > \
    --protein     < path to protein (aa) fasta file > \
    --ill         < path to directory that immediately contains all R1/R2 fastqs > \
    --iso         < path to directory that immediately contains all long read fastqs > \
    --masked      < bool > \
    --skip_qc     < bool > \
    --skip_trim   < bool > \
    --nlrs        < bool > \
    --ipscan      < path to interproscan directory > \
    --genemark    < path to prepared genemark directory (see `optional files` section) > \
    --design      < tsv file with expected gff files and weights for mikado > \
    --scoring     < path to mikado scoring file (e.g., plant.yaml) > \
    --homology    < path to homology fasta file (e.g., uniprot for you phylum )> \
    -profile local,four \
    -resume
```

Note the profile here is set up for use on a local server, but will likely require modification for your job. The `four` local profile is set up to use approximately 4 cpus maximum. Other presets exist in `conf/local.conf` and you can create your own by copying those examples.


## running on a slurm server

First, do you have a working nextflow/apptainer version?

```
nextflow -version
apptainer --version
```

Below is a sample sbatch script to run the pipeline. You'll need to replace the `<>` values with those that make sense for your use case.

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

nextflow  ~/nextflow/ragnarok/main.nf \
    --publish_dir < path to results location > \
    --genome      < path to genome in fasta > \
    --cds         < path to cds fasta file > \
    --protein     < path to protein (aa) fasta file > \
    --ill         < path to directory that immediately contains all R1/R2 fastqs > \
    --iso         < path to directory that immediately contains all long read fastqs > \
    --masked      < bool > \
    --skip_qc     < bool > \
    --skip_trim   < bool > \
    --nlrs        < bool > \
    --ipscan      < path to interproscan directory > \
    --genemark    < path to prepared genemark directory (see `optional files` section) > \
    --design      < tsv file with expected gff files and weights for mikado > \
    --scoring     < path to mikado scoring file (e.g., plant.yaml) > \
    --homology    < path to homology fasta file (e.g., uniprot for you phylum )> \
    -profile slurm,custom \
    -resume
```

Based on your qos/partition, you may want to modify the `conf/slurm.config` and `conf/slurm_custom.config` files to handle your dataset.

### a few notes on slurm qos and partitions

This pipeline currently relies on four qos/partition configurations on the UTK ISAAC-NG system.

- short : maximum 3 hours, 12 jobs submitted
- campus : maximum 1 day, 94 jobs submitted
- long : maximum 6 days, 14 jobs submitted
- gpu : maximum 4 hours, 1 job submitted

For each of the three labels found in `slurm.config` maxForks can be adjusted to qos that work on other sytems, and the imported clusterOptions for each can be updated with specific qos/partition/account information found in `slurm_custom.config`.

To check the limits for a given qos, replace short with any qos:

```
sacctmgr show qos where name=short
```

...or

```
scontrol show partition short
```

# development planning:

## desired input considerations

## necessary steps  
- [x] alignment 
  - [x] alignment - STAR
  - [x] alignment - long reads
    - [ ] consider options for pacbio (splice) v. nanopore
- [x] helixer
- [x] stringtie
- [x] gffread
- [x] miniprot
- [x] mikado2
- [x] busco
- [x] handle additional gff input paths
- [x] determine steps where copying to publish_dir is needed

## optional steps  
- [x] QC
- [x] trimming
- [x] FindPlantNLRs annotation
    - [x] FindPlantNLRs detection for parsing mikado2 weights
    - [x] FindPlantNLRs add to all_gffs_ch
- [x] allow non land_plant model for Helixer (opened issue)
- [ ] mikado2 plants.yaml max intron length (10k is good)
    - [ ] assets dir with angiosperms/gymnosperm preset
    - [ ] have a table of stats (toss citations in there)
- [x] allow for stringtie --mixed
- [x] allow "ill", "iso", or "mixed"
- [x] allow STAR and minimap to take in multiple fastq files

## obstacles/consider
- [x] mikado2 quay container is broken
- [ ] edta run fails on citrus genome (send error code) (rename scafs)
    - may need to rename scaffolds after running
- [ ] find braker3 logs for duplicated genes
    - full table tsv, busco id numbers, status
    - breakdown of types (gff files)
- [x] findplantnlrs, if interrupted, requires removal of FindPlantNLRs dir
```{bash}
# behavior fixed, but keep in case...
find work/ -type d -name "FindPlantNLRs"
```

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

# tools used in ragnarok

- [BUSCO](https://busco.ezlab.org)
- [compleasm](https://github.com/huangnengCSU/compleasm)
- [diamond](https://github.com/bbuchfink/diamond)
- [EDTA](https://github.com/oushujun/EDTA)
- [EnTAP2](https://entap.readthedocs.io/en/latest/index.html)
- [fastp](https://github.com/OpenGene/fastp)
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [FindPlantNLRs](https://github.com/ZhenyanLuo/FindPlantNLRs/tree/docker_version)
- [gffread](https://github.com/gpertea/gffread)
- [Helixer](https://github.com/weberlab-hhu/Helixer)
- [Mikado](https://mikado.readthedocs.io/en/stable/)
- [minimap2](https://github.com/lh3/minimap2)
- [miniprot](https://github.com/lh3/miniprot)
- [multiqc](https://github.com/MultiQC/MultiQC)
- [pandas](https://pandas.pydata.org)
- [samtools](https://www.htslib.org)
- [STAR](https://github.com/alexdobin/STAR)
- [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml)
- [TransDecoder](https://github.com/TransDecoder/TransDecoder)

## tool images

- busco:quay.io/biocontainers/busco:5.8.2--pyhdfd78af_0
- compleasm:quay.io/biocontainers/compleasm:0.2.6--pyh7cba7a3_0
- diamond:quay.io/biocontainers/diamond:2.1.11--h5ca1c30_1
- edta:quay.io/biocontainers/edta:2.2.2--hdfd78af_1
- entap:docker://plantgenomics/entap:2.2.0
- fastp:quay.io/biocontainers/fastp:0.23.4--h125f33a_5
- fastqc:quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
- findplantnlrs:docker://ryandk/findplantnlrs:latest
- gffread:quay.io/biocontainers/gffread:0.12.7--h077b44d_6
- helixer:docker://gglyptodon/helixer-docker:helixer_v0.3.4_cuda_12.2.2-cudnn8
- mikado2:docker://gemygk/mikado:v2.3.5rc2
- minimap2:quay.io/biocontainers/minimap2:2.28--h577a1d6_4
- miniprot:quay.io/biocontainers/miniprot:0.13--h577a1d6_2
- multiqc:quay.io/biocontainers/multiqc:1.24.1--pyhdfd78af_0
- pandas:quay.io/biocontainers/pandas:1.5.2
- samtools:quay.io/biocontainers/samtools:1.20--h50ea8bc_1
- star:quay.io/biocontainers/star:2.7.11a--h0033a41_0
- stringtie:quay.io/biocontainers/stringtie:3.0.0--h29c0135_0
- transdecoder:quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0

See conf/containers.config for most current versions.
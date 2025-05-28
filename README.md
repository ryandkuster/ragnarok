# Ragnarok

```
################################################################################
#                                   RAGNAROK!                                  #
#                         RApid Genome anNotAtion ROcKs!                       #
#                          Chris Gottschalk 12/6/2024                          #
#                                    v2.1                                      #
################################################################################
```

# requirements

## software
nextflow (22.10.4+)  
apptainer (1.1.8+)

## required files

|parameter|type|description|
|:-|:-|:-|
|`--genome`|.fna|Provide your assembly (ideally with *simple* names if using EDTA masking option).|
|`--ill`|string|Required if `iso` not used. Path to a directory containing paired end fastq files. Can end with directory name (Ex: path/to/files/) or a specific prefix of the paired files (Ex: path/to/files/reads_P1). Using specific prefix name will search for the paired end read files ending in R1.[fq|fastq].gz .|
|`--iso`|string|Required if `ill` not used. Path to a directory containing long read fastq files. Can end with directory name (Ex: path/to/files/) or a specific prefix of the paired files (Ex: path/to/files/reads_P1). Using specific prefix name will search for the paired end read files ending in [fq|fastq].gz .|
|`--protein`|.faa|Protein file for miniprot (e.g., closest ref species).|
|`--homology`|.faa|Mikado protein homology file (e.g., [uniprot 33090 for viridiplantae](https://www.uniprot.org/uniprotkb?query=viridiplantae&facets=reviewed%3Atrue)).|
|`--scoring`|.yaml|Mikado scoring file (e.g., [plant.yaml](https://github.com/EI-CoreBioinformatics/mikado/tree/master/Mikado/configuration/scoring_files)).|
|`--design`|.tsv|Mikado configuration table file (see Mikado [documentation](https://mikado.readthedocs.io/en/stable/Tutorial/) and below).|

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
> _The first field is intentionally missing as Ragnarok will produce these outputs._

Example tsv configuration for `--nlrs true` (assets/mikado_nlr_conf.tsv):

```
	hx	True		False	False
	st	True	1	False	True
	tr	False	-0.5	False	False
	mp	True	1	False	False
	nlr	True	1	False	False
```

> [!NOTE]
> _The first field is intentionally missing as Ragnarok will produce these outputs._

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
> _The filepath to existing gffs will need to be provided._


## optional files

|parameter|type|description|
|:-|:-|:-|
|`--cds`|.fna|CDS file for your species for use with `--perform_masking true` (used by EDTA)|
|`--ipscan`|directory|Locally stored interproscan for use with `--nlrs true` ([64-bit download](https://www.ebi.ac.uk/interpro/download/InterProScan/))|
|`--genemark`|directory|Genemark with key configured for use with `--nlrs true` (see assets/genemark_setup.sh)|

## additional parameters

|parameter|type|description|default|
|:-|:-|:-|:-|
|`--lineage`|url|URL path to helixer model https://zenodo.org/records/10836346.|"https://zenodo.org/records/10836346/files/land_plant_v0.3_a_0080.h5"|
|`--subseq_len`|int|Helixer subseq length.|64152|
|`--skip_qc`|bool|Perform fastqc/multiqc on raw read data.|true|
|`--skip_trim`|bool|Perform adapter trimming on raw read data|true|
|`--minimum_length`|int|Use with `--skip_trim`, minimum length read to keep when adapter trimming.|50|
|`--perform_masking`|bool|Run EDTA to mask input genome (recommended).|false|
|`--masking_threshold`|int|Use with `perform_masking` to custom hard-mask TEanno models >= this length.|1000|
|`--skip_st`|bool|Requires st, tr, mp gff files locally (in `--design` file), bypass Stringtie steps.|false|
|`--skip_hx`|bool|Requires hx file locally (in `--design` file), bypass Helixer step.|false|
|`--contam`|comma-sep list|Entap list of taxa to consider as contaminants.|"insecta,fungi,bacteria"|

# getting started

## running on a local server

First, do you have a working nextflow/apptainer version?

```
nextflow -version
apptainer --version
```

Below is a sample script to run the pipeline. You'll need to replace the `<>` values with those that make sense for your use case.

```
nextflow run ~/nextflow/ragnarok/main.nf \
    --publish_dir     < path to results location > \
    --genome          < path to genome in fasta > \
    --cds             < path to cds fasta file > \
    --protein         < path to protein (aa) fasta file > \
    --ill             < path to directory that immediately contains all R1/R2 fastqs > \
    --iso             < path to directory that immediately contains all long read fastqs > \
    --perform_masking < bool > \
    --skip_qc         < bool > \
    --skip_trim       < bool > \
    --nlrs            < bool > \
    --ipscan          < path to interproscan directory > \
    --genemark        < path to prepared genemark directory (see `optional files` section) > \
    --design          < tsv file with expected gff files and weights for mikado > \
    --scoring         < path to mikado scoring file (e.g., plant.yaml) > \
    --homology        < path to homology fasta file (e.g., uniprot for you phylum )> \
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

nextflow run ~/nextflow/ragnarok/main.nf \
    --publish_dir     < path to results location > \
    --genome          < path to genome in fasta > \
    --cds             < path to cds fasta file > \
    --protein         < path to protein (aa) fasta file > \
    --ill             < path to directory that immediately contains all R1/R2 fastqs > \
    --iso             < path to directory that immediately contains all long read fastqs > \
    --perform_masking < bool > \
    --skip_qc         < bool > \
    --skip_trim       < bool > \
    --nlrs            < bool > \
    --ipscan          < path to interproscan directory > \
    --genemark        < path to prepared genemark directory (see `optional files` section) > \
    --design          < tsv file with expected gff files and weights for mikado > \
    --scoring         < path to mikado scoring file (e.g., plant.yaml) > \
    --homology        < path to homology fasta file (e.g., uniprot for you phylum )> \
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
- [ ] QC
    - [x] FastQC for ill input
    - [ ] LongQC for iso input
- [ ] trimming
    - [x] fastp for ill input
    - [ ] ? for iso input
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
- [ ] add EnTAP2 functional annotations
- [ ] add liftover input (gff + genome) from closely related accession

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
    v40["lineage"]
    v56["scoring"]
    v57["homology"]
    v2["iso"]
    v53["genemark"]
    v0["ill"]
    v51["ipscan"]
    v46["nlrs"]
    v4["skip_qc"]
    v8["minimum_length"]
    v39["skip_hx"]
    v14["genome"]
    v15["cds"]
    v34["protein"]
    v42["subseq_len"]
    v45["design"]
    v7["skip_trim"]
    v13["perform_masking"]
    v25["skip_st"]
    end
    v5([FASTQC_RAW])
    v6([MULTIQC_RAW])
    v9([FASTP_ADAPTERS])
    v11([FASTQC_TRIM])
    v12([MULTIQC_TRIM])
    v16([EDTA])
    v19([STAR_INDEX_NA])
    v20([STAR_MAP])
    v22([SAM_SORT])
    v23([MINIMAP2])
    v24([SAM_SORT_LONG])
    v26([STRINGTIE_MIX])
    v28([STRINGTIE])
    v30([STRINGTIE])
    v32([GFFREAD])
    v33([TRANSDECODER])
    v35([MINIPROT])
    v41([HELIXER_DB])
    v43([HELIXER])
    v47([PARSE_INPUT])
    v50([FPNLRS_SETUP])
    v52([FINDPLANTNLRS])
    v54([ANNOTATENLRS])
    v58([MIKADO_CONF])
    v59([TRANSDECODER_ORF])
    v60([DIAMOND])
    v61([THE_GRANDMASTER])
    v62([GFFREAD_FINAL])
    v63([BUSCO])
    v64([COMPLEASM_DB])
    v65([COMPLEASM])
    v0 --> v5
    v5 --> v6
    v0 --> v9
    v8 --> v9
    v9 --> v11
    v11 --> v12
    v14 --> v16
    v15 --> v16
    v14 --> v19
    v19 --> v20
    v9 --> v20
    v20 --> v22
    v2 --> v23
    v14 --> v23
    v23 --> v24
    v22 --> v26
    v24 --> v26
    v22 --> v28
    v24 --> v30
    v14 --> v32
    v30 --> v32
    v14 --> v33
    v30 --> v33
    v33 --> v35
    v34 --> v35
    v14 --> v35
    v40 --> v41
    v41 --> v43
    v42 --> v43
    v14 --> v43
    v39 --> v47
    v25 --> v47
    v45 --> v47
    v46 --> v47
    v14 --> v50
    v50 --> v52
    v51 --> v52
    v51 --> v54
    v52 --> v54
    v53 --> v54
    v32 --> v58
    v33 --> v58
    v35 --> v58
    v54 --> v58
    v56 --> v58
    v57 --> v58
    v43 --> v58
    v14 --> v58
    v47 --> v58
    v58 --> v59
    v57 --> v60
    v58 --> v60
    v57 --> v61
    v58 --> v61
    v59 --> v61
    v60 --> v61
    v14 --> v61
    v61 --> v62
    v14 --> v62
    v62 --> v63
    v64 --> v65
    v62 --> v65
    end
```

# tools used in ragnarok

- [bedtools](https://bedtools.readthedocs.io/en/latest/)
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

- agat:quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0
- bedtools:quay.io/biocontainers/bedtools:2.31.1--h13024bc_3
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

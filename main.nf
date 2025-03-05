nextflow.enable.dsl=2

// adapter removal
include { FASTP_ADAPTERS          } from "./modules/fastp.nf"

// qc
include { FASTQC as FASTQC_RAW    } from "./modules/fastqc.nf"
include { FASTQC as FASTQC_TRIM   } from "./modules/fastqc.nf"
include { MULTIQC as MULTIQC_RAW  } from "./modules/multiqc.nf"
include { MULTIQC as MULTIQC_TRIM } from "./modules/multiqc.nf"

// star + samtools
include { STAR_INDEX_NA           } from "./modules/star.nf"
include { STAR_MAP                } from "./modules/star.nf"
include { SAM_SORT                } from "./modules/samtools.nf"
include { SAM_INDEX               } from "./modules/samtools.nf"

// others
include { EDTA                    } from "./modules/edta.nf"
include { STRINGTIE               } from "./modules/stringtie.nf"
include { GFFREAD                 } from "./modules/gffread.nf"
include { TRANSDECODER            } from "./modules/transdecoder.nf"
include { MINIPROT                } from "./modules/miniprot.nf"
include { HELIXER                 } from "./modules/helixer.nf"
include { MIKADO2                 } from "./modules/mikado2.nf"

workflow {

    /*
    --------------------------------------------------------------------
        parse inputs and prepare indices and interval files
    --------------------------------------------------------------------
    */

    fastq_ch = Channel.fromFilePairs("${params.fastq_pe}/*_R{1,2}*", checkIfExists: true, flat:true)

    /*
    --------------------------------------------------------------------
        adapter removal and qc summary stats
    --------------------------------------------------------------------
    */

    if (!params.skip_qc) {
        FASTQC_RAW(
            fastq_ch,
            "raw")
        MULTIQC_RAW(
            FASTQC_RAW.out.fastq_ch.collect(),
            "raw")
    }

    if (!params.skip_trim) {
        FASTP_ADAPTERS(fastq_ch, params.minimum_length)
        fastq_ch = FASTP_ADAPTERS.out.reads

        if (!params.skip_qc) {
            FASTQC_TRIM(fastq_ch, "trimmed")
            MULTIQC_TRIM(FASTQC_TRIM.out.fastq_ch.collect(), "trimmed")
        }
    }

    /*
    --------------------------------------------------------------------
        read alignment
    --------------------------------------------------------------------
    */

    if (params.masked == false) {
        EDTA(params.genome,
             params.cds)
    }

    if (params.read_type == "ill") {
 
    STAR_INDEX_NA(params.genome)
 
    STAR_MAP(fastq_ch,
             STAR_INDEX_NA.out.star_idx)
    }
 
    SAM_SORT(STAR_MAP.out.bam_ch)

    // TODO: consider allowing bam as input to bypass star mapping
    STRINGTIE(SAM_SORT.out.sort_ch)
 
    GFFREAD(STRINGTIE.out.gtf_ch,
            params.genome)
 
    TRANSDECODER(STRINGTIE.out.gtf_ch,
            params.genome)
 
    MINIPROT(TRANSDECODER.out.est_ch,
            params.genome,
            params.protein)

    HELIXER(params.genome,
            params.lineage)

    /*
    --------------------------------------------------------------------
        mikado2
    --------------------------------------------------------------------
    */
}

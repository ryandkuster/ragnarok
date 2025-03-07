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
include { HELIXER_DB              } from "./modules/helixer.nf"
include { PARSE_INPUT             } from './modules/parse_input.nf'
include { MIKADO_CONF             } from "./modules/mikado2.nf"
include { DIAMOND                 } from './modules/diamond.nf'
include { THE_GRANDMASTER         } from './modules/mikado2.nf'
include { TRANSDECODER_ORF        } from './modules/transdecoder.nf'
include { GFFREAD_FINAL           } from './modules/gffread.nf'
include { BUSCO                   } from './modules/busco.nf'
include { COMPLEASM } from './modules/compleasm.nf'

workflow {

    /*
    --------------------------------------------------------------------
        get fastq inputs
    --------------------------------------------------------------------
    */

    fastq_ch = Channel.fromFilePairs("${params.fastq_pe}/*_R{1,2}*", checkIfExists: true, flat:true)

    /*
    --------------------------------------------------------------------
        adapter removal and qc summary stats (optional)
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
        SAM_SORT(STAR_MAP.out.bam_ch)
    }
 
    if (params.read_type == "iso") {}

    // TODO: consider allowing bam as input to bypass star mapping
    if (!params.skip_st) {
        STRINGTIE(SAM_SORT.out.sort_ch)
 
        GFFREAD(STRINGTIE.out.st_ch,
                params.genome)
 
        TRANSDECODER(STRINGTIE.out.st_ch,
                     params.genome)
 
        MINIPROT(TRANSDECODER.out.tr_ch,
                 params.genome,
                 params.protein)

        HELIXER_DB(params.lineage)

        st_gff = GFFREAD.out.st_gff_ch
        tr_gff = TRANSDECODER.out.tr_ch
        mp_gff = MINIPROT.out.mp_ch
    }

    if (!params.skip_hx) {
        HELIXER(params.genome,
                HELIXER_DB.out.db_ch,
                params.subseq_len)

        hx_gff = HELIXER.out.hx_ch
    }

    /*
    --------------------------------------------------------------------
        mikado2
    --------------------------------------------------------------------
    */

    PARSE_INPUT(params.design,
                params.skip_st,
                params.skip_hx)

    PARSE_INPUT.out.path_ch
        .splitCsv( header: false, sep: ',' )
        .set { gff_path_ch }

    gff_path_ch.concat(hx_gff, st_gff, tr_gff, mp_gff).set{all_gff_ch}
    all_gff_ch.view()

    MIKADO_CONF(all_gff_ch.collect(),
                    PARSE_INPUT.out.design_ch,
                    params.genome,
                    params.scoring,
                    params.homology)

    TRANSDECODER_ORF(MIKADO_CONF.out.mk_ch)

    DIAMOND(params.homology,
            MIKADO_CONF.out.mk_ch)

    THE_GRANDMASTER(MIKADO_CONF.out.yaml_ch,
                    MIKADO_CONF.out.mk_ch,
                    DIAMOND.out.dmnd_ch,
                    TRANSDECODER_ORF.out.orf_ch,
                    params.homology,
                    params.genome)

    GFFREAD_FINAL(THE_GRANDMASTER.out.gm_ch,
                  params.genome)

    BUSCO(GFFREAD_FINAL.out.final_ch)
    COMPLEASM(GFFREAD_FINAL.out.final_ch)
}

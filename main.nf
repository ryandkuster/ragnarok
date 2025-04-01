nextflow.enable.dsl=2

// adapter removal
include { FASTP_ADAPTERS            } from "./modules/fastp.nf"

// qc
include { FASTQC as FASTQC_RAW      } from "./modules/fastqc.nf"
include { FASTQC as FASTQC_TRIM     } from "./modules/fastqc.nf"
include { MULTIQC as MULTIQC_RAW    } from "./modules/multiqc.nf"
include { MULTIQC as MULTIQC_TRIM   } from "./modules/multiqc.nf"

// star + samtools
include { STAR_INDEX_NA             } from "./modules/star.nf"
include { STAR_MAP                  } from "./modules/star.nf"
include { SAM_SORT                  } from "./modules/samtools.nf"
include { SAM_SORT as SAM_SORT_LONG } from "./modules/samtools.nf"
include { SAM_INDEX                 } from "./modules/samtools.nf"

// others
include { EDTA                      } from "./modules/edta.nf"
include { STRINGTIE                 } from "./modules/stringtie.nf"
include { STRINGTIE_MIX             } from "./modules/stringtie.nf"
include { GFFREAD                   } from "./modules/gffread.nf"
include { TRANSDECODER              } from "./modules/transdecoder.nf"
include { MINIPROT                  } from "./modules/miniprot.nf"
include { HELIXER                   } from "./modules/helixer.nf"
include { HELIXER_DB                } from "./modules/helixer.nf"
include { PARSE_INPUT               } from './modules/parse_input.nf'
include { MIKADO_CONF               } from "./modules/mikado2.nf"
include { DIAMOND                   } from './modules/diamond.nf'
include { THE_GRANDMASTER           } from './modules/mikado2.nf'
include { TRANSDECODER_ORF          } from './modules/transdecoder.nf'
include { GFFREAD_FINAL             } from './modules/gffread.nf'
include { BUSCO                     } from './modules/busco.nf'
include { COMPLEASM                 } from './modules/compleasm.nf'
include { COMPLEASM_DB              } from './modules/compleasm.nf'
include { FPNLRS_SETUP              } from './modules/findplantnlrs.nf'
include { FINDPLANTNLRS             } from './modules/findplantnlrs.nf'
include { ANNOTATENLRS              } from './modules/findplantnlrs.nf'
include { MINIMAP2                  } from './modules/minimap2.nf'

workflow {

    /*
    --------------------------------------------------------------------
        get fastq inputs
    --------------------------------------------------------------------
    */

    if (params.ill) {
        fastq_ch = Channel.fromFilePairs("${params.ill}/*_{R,}{1,2}*", checkIfExists: true, flat:true)
    }
    if (params.iso) {
        long_ch = Channel.fromPath("${params.iso}/*{fq,fastq}*", checkIfExists: true)
        long_ch.view()
    }

    /*
    --------------------------------------------------------------------
        adapter removal and qc summary stats (optional)
    --------------------------------------------------------------------
    */

    if (!params.skip_qc) {
        if (params.ill) {
            FASTQC_RAW(
                fastq_ch,
                "raw")
            MULTIQC_RAW(
                FASTQC_RAW.out.fastq_ch.collect(),
                "raw")
        }
    }

    // TODO: add QC for iso (long_ch) reads
    if (!params.skip_trim) {
        if (params.ill) {
            FASTP_ADAPTERS(fastq_ch, params.minimum_length)
            fastq_ch = FASTP_ADAPTERS.out.reads

            if (!params.skip_qc) {
                FASTQC_TRIM(fastq_ch, "trimmed")
                MULTIQC_TRIM(FASTQC_TRIM.out.fastq_ch.collect(), "trimmed")
            }
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

    if (params.ill) {
        STAR_INDEX_NA(params.genome)
        STAR_MAP(
            fastq_ch,
            STAR_INDEX_NA.out.star_idx)
        STAR_MAP.out.bam_ch.map { item -> [item[1]] }.collect().set { star_ch }
        SAM_SORT(
            star_ch.collect(),
            "short")
    }
 
    if (params.iso) {
        MINIMAP2(
            params.genome,
            long_ch.collect())
        SAM_SORT_LONG(
            MINIMAP2.out.mp_ch,
            "long")
    }

    // TODO: consider allowing bam as input to bypass star mapping
    if (!params.skip_st) {

        if (params.ill && params.iso){
            STRINGTIE_MIX(
                SAM_SORT.out.sort_ch,
                SAM_SORT_LONG.out.sort_ch)
            STRINGTIE_MIX.out.st_ch.set{ st_ch }
        } else if (params.ill) {
            STRINGTIE(SAM_SORT.out.sort_ch)
            STRINGTIE.out.st_ch.set{ st_ch }
        } else if (params.iso) {
            STRINGTIE(SAM_SORT_LONG.out.sort_ch)
            STRINGTIE.out.st_ch.set{ st_ch }
        }
 
        GFFREAD(st_ch,
                params.genome)
 
        TRANSDECODER(st_ch,
                     params.genome)
 
        MINIPROT(TRANSDECODER.out.tr_ch,
                 params.genome,
                 params.protein)

        st_gff = GFFREAD.out.st_gff_ch
        tr_gff = TRANSDECODER.out.tr_ch
        mp_gff = MINIPROT.out.mp_ch
    }

    if (!params.skip_hx) {
        HELIXER_DB(params.lineage)
        HELIXER(params.genome,
                HELIXER_DB.out.db_ch,
                params.subseq_len)

        hx_gff = HELIXER.out.hx_ch
    }

    /*
    --------------------------------------------------------------------
        parse user design and gather input gffs
    --------------------------------------------------------------------
    */

    PARSE_INPUT(params.design,
                params.skip_st,
                params.skip_hx,
                params.nlrs)

    PARSE_INPUT.out.path_ch
        .splitCsv( header: false, sep: ',' )
        .set { gff_path_ch }

    gff_path_ch.concat(hx_gff, st_gff, tr_gff, mp_gff).set{ all_gff_ch }

    /*
    --------------------------------------------------------------------
        find plant NLRs
    --------------------------------------------------------------------
    */

    if (params.nlrs == true) {
        FPNLRS_SETUP(params.genome)
        FINDPLANTNLRS(FPNLRS_SETUP.out.fpnlr_db_ch,
                      params.ipscan)
        ANNOTATENLRS(FINDPLANTNLRS.out.fpnlr_ch,
                      params.ipscan,
                      params.genemark)
        all_gff_ch.mix(ANNOTATENLRS.out.nlr_gff_ch).set{ all_gff_ch }
    }

    /*
    --------------------------------------------------------------------
        mikado2
    --------------------------------------------------------------------
    */

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

    /*
    --------------------------------------------------------------------
        BUSCO quality metrics
    --------------------------------------------------------------------
    */

    BUSCO(GFFREAD_FINAL.out.final_ch)
    COMPLEASM_DB()
    COMPLEASM(GFFREAD_FINAL.out.final_ch,
              COMPLEASM_DB.out.db_ch)

}

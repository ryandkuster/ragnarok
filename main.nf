nextflow.enable.dsl=2

// raw read qc and adapter removal
include { FASTQC as FASTQC_RAW      } from "./modules/fastqc.nf"
include { FASTQC as FASTQC_TRIM     } from "./modules/fastqc.nf"
include { MULTIQC as MULTIQC_RAW    } from "./modules/multiqc.nf"
include { MULTIQC as MULTIQC_TRIM   } from "./modules/multiqc.nf"
include { FASTP_ADAPTERS            } from "./modules/fastp.nf"

// alignment
include { STAR_INDEX_NA             } from "./modules/star.nf"
include { STAR_MAP                  } from "./modules/star.nf"
include { SAM_SORT                  } from "./modules/samtools.nf"
include { SAM_SORT as SAM_SORT_LONG } from "./modules/samtools.nf"
include { SAM_INDEX                 } from "./modules/samtools.nf"
include { MINIMAP2                  } from './modules/minimap2.nf'

// masking
include { EDTA                      } from "./modules/edta.nf"
include { EDTA_THRESHOLD            } from "./modules/edta.nf"

// main annotations
include { STRINGTIE                 } from "./modules/stringtie.nf"
include { STRINGTIE_MIX             } from "./modules/stringtie.nf"
include { GFFREAD                   } from "./modules/gffread.nf"
include { TRANSDECODER              } from "./modules/transdecoder.nf"
include { MINIPROT                  } from "./modules/miniprot.nf"
include { HELIXER                   } from "./modules/helixer.nf"
include { HELIXER_DB                } from "./modules/helixer.nf"

// intermediates
include { PARSE_INPUT               } from './modules/parse_input.nf'
include { SCAF2NUM                  } from './modules/parse_input.nf'
include { NUM2SCAF as NUM2SCAF_1    } from './modules/parse_input.nf'
include { NUM2SCAF as NUM2SCAF_2    } from './modules/parse_input.nf'
include { NOZIP_REF                 } from './modules/parse_input.nf'
include { PROT_FIX                  } from './modules/parse_input.nf'
include { MIKADO_CONF               } from "./modules/mikado2.nf"
include { DIAMOND                   } from './modules/diamond.nf'
include { THE_GRANDMASTER           } from './modules/mikado2.nf'
include { TRANSDECODER_ORF          } from './modules/transdecoder.nf'
include { GFFREAD_MIKADO            } from './modules/gffread.nf'
include { GFFREAD_ENTAP             } from './modules/gffread.nf'
include { AGAT_SUBSET               } from './modules/agat.nf'

// annotation qc
include { BUSCO                     } from './modules/busco.nf'
include { COMPLEASM                 } from './modules/compleasm.nf'
include { COMPLEASM_DB              } from './modules/compleasm.nf'

// find plant nlrs
include { LIFTOFF                   } from './modules/liftover.nf'

// find plant nlrs
include { FPNLRS_SETUP              } from './modules/findplantnlrs.nf'
include { FINDPLANTNLRS             } from './modules/findplantnlrs.nf'
include { ANNOTATENLRS              } from './modules/findplantnlrs.nf'

// find plant nlrs
include { ENTAP_INI                 } from './modules/entap.nf'
include { ENTAP_RUN                 } from './modules/entap.nf'

workflow {

    /*
    --------------------------------------------------------------------
        parse user design and gather input gffs
    --------------------------------------------------------------------
    */

    PARSE_INPUT(params.design,
                params.skip_st,
                params.skip_hx,
                params.nlrs,
                params.lo_genome,
                params.lo_gff)

    PARSE_INPUT.out.path_ch
        .splitCsv( header: false, sep: ',' )
        .set { gff_path_ch }

    NOZIP_REF(params.genome)
    NOZIP_REF.out.nozip_ref_ch.set{ mk_genome }

    /*
    --------------------------------------------------------------------
        get fastq inputs
    --------------------------------------------------------------------
    */

    if (params.ill) {
        fastq_ch = Channel.fromFilePairs("${params.ill}*_{R,}{1,2}*", checkIfExists: true, flat:true)
    }
    if (params.iso) {
        long_ch = Channel.fromPath("${params.iso}*{fq,fastq}*", checkIfExists: true)
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
        genome masking
    --------------------------------------------------------------------
    */

    if (params.perform_masking == true) {
        // Rename the input fastas with simple naming scheme.
        SCAF2NUM(
            params.genome,
            params.cds)

        // Perform masking of repetitive regions.
        EDTA(
            SCAF2NUM.out.s2n_ref_ch,
            SCAF2NUM.out.s2n_cds_ch)

        // Convert the EDTA hard-masked fasta to original naming scheme.
        NUM2SCAF_1(
            EDTA.out.mask_ch,
            SCAF2NUM.out.s2n_json_ch)

        // Use the EDTA hard-masked genome for stringtie steps.
        NUM2SCAF_1.out.n2s_ref_ch.set{ st_genome }

        // Optionally, allow a custom TE anno size (default 1kb) to
        // allow more or less stringency for helixer annotations.
        if (params.masking_threshold != false) {
            EDTA_THRESHOLD(
                EDTA.out.edta_ch,
                SCAF2NUM.out.s2n_ref_ch,
                params.masking_threshold)
            NUM2SCAF_2(
                EDTA_THRESHOLD.out.threshold_ch,
                SCAF2NUM.out.s2n_json_ch)
            NUM2SCAF_2.out.n2s_ref_ch.set{ hx_genome }
        } else {
            NUM2SCAF_1.out.n2s_ref_ch.set{ hx_genome }
        }

    } else {
        st_genome = Channel.fromPath("${params.genome}", checkIfExists: true)
        hx_genome = Channel.fromPath("${params.genome}", checkIfExists: true)
    }

    /*
    --------------------------------------------------------------------
        read alignment
    --------------------------------------------------------------------
    */

    if (params.ill) {
        STAR_INDEX_NA(st_genome)

        // star_idx must be a value channel
        star_idx = STAR_INDEX_NA.out.star_idx.collect()

        STAR_MAP(
            fastq_ch,
            star_idx)

        STAR_MAP.out.bam_ch.map { item -> [item[1]] }.collect().set { star_ch }
        SAM_SORT(
            star_ch.collect(),
            "short")
    }
 
    if (params.iso) {
        MINIMAP2(
            st_genome,
            long_ch.collect())
        SAM_SORT_LONG(
            MINIMAP2.out.mp_ch,
            "long")
    }

    /*
    --------------------------------------------------------------------
        core ragnarok annotation steps
    --------------------------------------------------------------------
    */

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
                st_genome)
 
        TRANSDECODER(st_ch,
                     st_genome)
 
        MINIPROT(TRANSDECODER.out.tr_ch,
                 st_genome,
                 params.protein)

        st_gff = GFFREAD.out.st_gff_ch
        tr_gff = TRANSDECODER.out.tr_ch
        mp_gff = MINIPROT.out.mp_ch
    }

    if (!params.skip_hx) {
        HELIXER_DB(params.lineage)
        HELIXER(hx_genome,
                HELIXER_DB.out.db_ch,
                params.subseq_len)

        hx_gff = HELIXER.out.hx_ch
    }

    if (params.skip_hx == false && params.skip_st == false) {
        gff_path_ch.concat(hx_gff, st_gff, tr_gff, mp_gff).set{ all_gff_ch }
    } else if (params.skip_hx == false && params.skip_st == true) {
        gff_path_ch.concat(hx_gff).set{ all_gff_ch }
    } else if (params.skip_hx == true && params.skip_st == false) {
        gff_path_ch.concat(st_gff, tr_gff, mp_gff).set{ all_gff_ch }
    } else if (params.skip_hx == true && params.skip_st == true) {
        gff_path_ch.set{ all_gff_ch }
    }

    /*
    --------------------------------------------------------------------
        lifton (optional)
    --------------------------------------------------------------------
    */

    if (params.lo_genome) {
        LIFTOFF(mk_genome,
               params.lo_genome,
               params.lo_gff)
        all_gff_ch.mix(LIFTOFF.out.lo_ch).set{ all_gff_ch }
    }

    /*
    --------------------------------------------------------------------
        find plant NLRs (optional)
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
                    mk_genome,
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
                    mk_genome)

    GFFREAD_MIKADO(THE_GRANDMASTER.out.gm_ch,
                  mk_genome)

    /*
    --------------------------------------------------------------------
        EnTAP functional annotation
    --------------------------------------------------------------------
    */

    entap_conf = file("${projectDir}/assets/template_entap_config.ini")
    entap_run = file("${projectDir}/assets/template_entap_run.params")
    ENTAP_INI(entap_conf,
              entap_run)
    PROT_FIX(GFFREAD_MIKADO.out.final_ch)
    ENTAP_RUN(ENTAP_INI.out.conf_ch,
              ENTAP_INI.out.run_ch,
              ENTAP_INI.out.db_ch,
              PROT_FIX.out.prot_ch)

    AGAT_SUBSET(THE_GRANDMASTER.out.gm_ch,
                ENTAP_RUN.out.annot_ch,
                params.final_prefix)

    /*
    --------------------------------------------------------------------
        BUSCO quality metrics
    --------------------------------------------------------------------
    */

    GFFREAD_ENTAP(AGAT_SUBSET.out.subset_pass_ch,
                  mk_genome,
                  params.final_prefix)
    BUSCO(GFFREAD_ENTAP.out.final_ch,
          params.busco_db,
          params.final_prefix)
    COMPLEASM_DB(params.busco_db)
    COMPLEASM(GFFREAD_ENTAP.out.final_ch,
              COMPLEASM_DB.out.db_ch,
              params.busco_db,
              params.final_prefix)

}

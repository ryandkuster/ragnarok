process MIKADO2 {
    label 'mikado2'
    label 'campus'

    time 10.h
    cpus 20
    memory 4.GB

    input:
        tuple val(sample), path(gff)
        path genome
        path protein

    output:
        tuple val(sample), path("*_aa_miniprot.gff"), emit: mk_ch

    script:
        """
        echo "M_ioensis_hap1_helixer_rd1.gff3   hx  True        False   False"
        echo "stringtie_LR.gtf  st  True    1   False   True"
        echo "transcripts.fasta.transdecoder.genome.gff3    tr  False   -0.5    False   False"
        echo "fusca_aa_miniprot.gff mp  True    1   False   False"
        """
}


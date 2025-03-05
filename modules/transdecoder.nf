process TRANSDECODER {
    label 'transdecoder'
    label 'campus'

    time 2.h
    cpus 4
    memory 4.GB

    input:
        tuple val(sample), path(gtf)
        path genome

    output:
        tuple val(sample), path("*.genome.gff3"), emit: est_ch

    script:
        """
        gtf_genome_to_cdna_fasta.pl $gtf $genome > transcripts.fasta
        TransDecoder.LongOrfs -t transcripts.fasta
        TransDecoder.Predict -t transcripts.fasta
        gtf_to_alignment_gff3.pl $gtf > transcripts.gff3

        cdna_alignment_orf_to_genome_orf.pl \
            transcripts.fasta.transdecoder.gff3 \
            transcripts.gff3 \
            transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
        """
}


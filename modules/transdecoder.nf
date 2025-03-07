process TRANSDECODER {
    label 'transdecoder'
    label 'campus'

    time 2.h
    cpus 4
    memory 4.GB

    input:
        path(gtf)
        path genome

    output:
        path("*.genome.gff3"), emit: tr_ch

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

process TRANSDECODER_ORF {
    label 'transdecoder'
    label 'short'

    time 3.h
    cpus 4
    memory 20.GB

    input:
        tuple path(mk_gtf), path(mk_fasta)

    output:
        path("mikado_prepared.fasta.transdecoder.gff3"), emit: orf_ch

    script:
        """
        TransDecoder.LongOrfs -t $mk_fasta
        mv mikado_prepared.fasta.transdecoder_dir/longest_orfs.gff3 ./mikado_prepared.fasta.transdecoder.gff3
        """
}
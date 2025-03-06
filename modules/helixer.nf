process HELIXER_DB {
    label 'helixer'
    label 'short'

    time 1.h
    cpus 4
    memory 4.GB

    input:
        val(lineage)

    output:
        path("*.h5"), emit: db_ch

    script:
        """
        wget $lineage
        """
}

process HELIXER {
    label 'helixer'
    label 'gpu'

    time 4.h
    cpus 12
    memory 10.GB

    input:
        path(genome)
        path(model)
        val(subseq_len)

    output:
        path("*.gff3"), emit: hx_ch

    script:
        """
        Helixer.py \
            --model-filepath $model \
            --subsequence-length $subseq_len \
            --fasta-path $genome \
            --species helixer_species \
            --gff-output-path helixer.gff3
        """
}
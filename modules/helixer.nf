process HELIXER {
    label 'helixer'
    label 'gpu'

    time 24.h
    cpus 24
    memory 50.GB

    input:
        path(genome)
        val(lineage)

    output:
        path("*.gff3"), emit: gff_ch

    script:
        """
        Helixer.py \
            --lineage $lineage \
            --fasta-path $genome \
            --species helixer_species \
            --gff-output-path helixer.gff3
        """
}


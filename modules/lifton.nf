process LIFTON {
    label 'lifton'
    label 'campus'

    time 12.h
    cpus 12
    memory 10.GB

    input:
        path(genome)
        path(lo_genome)
        path(lo_gff)

    output:
        path("lifton.gff3"), emit: lo_ch

    script:
        """
        lifton \
            -g $lo_gff \
            -o lifton.gff3 \
            -copies \
            -sc 0.95 \
            $genome \
            $lo_genome
        """
}
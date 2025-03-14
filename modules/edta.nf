process EDTA {
    label 'edta'
    label 'long'

    time 48.h
    cpus 24
    memory 40.GB

    input:
        path(genome)
        path(cds)

    script:
        """
        # names have to be short!!! (<=13 chars)
        EDTA.pl \
            --genome $genome \
            --cds $cds \
            --overwrite 1 \
            --sensitive 1 \
            --anno 1 \
            --force 1 \
            --threads ${task.cpus}
        """
}


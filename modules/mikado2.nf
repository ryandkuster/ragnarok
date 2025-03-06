process THE_GRANDMASTER {
    label 'mikado2'
    label 'campus'

    time 10.h
    cpus 20
    memory 4.GB

    input:
        path(all_gffs)
        path(design_ch)
        path(genome)
        path(scoring)
        path(homology)

    script:
        """
        mikado configure \
            --list mikado.tsv \
            --reference $genome \
            --mode permissive \
            --scoring $scoring \
            -bt $homology \
            configuration.yaml
        mikado prepare --json-conf configuration.yaml
        """
}


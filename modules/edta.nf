process EDTA {
    label 'edta'
    label 'long'

    time 48.h
    cpus 24
    memory 100.GB

    publishDir(path: "${publish_dir}/edta_masking", mode: "copy")

    input:
        path(genome)
        path(cds)

    output:
        path("*.{masked,fa,gff3,pdf,sum}"), emit: edta_ch
        path("*.masked"), emit: mask_ch

    script:
        """
        # names have to be short!!! (<=13 chars)

        unset -f which

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

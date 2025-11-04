process EDTA {
    label 'edta'
    label 'long'

    time 48.h
    cpus 24
    memory 100.GB

    publishDir(path: "${params.publish_dir}/publish/edta_masking", mode: "copy")

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

process EDTA_THRESHOLD {
    label 'bedtools'
    label 'short'

    time 3.h
    cpus 24
    memory 10.GB

    publishDir(path: "${params.publish_dir}/publish/edta_masking", mode: "copy")

    input:
        path("*")
        path(genome)
        val(threshold)

    output:
        path("custom.masked.fna"), emit: threshold_ch

    script:
        """
        grep -v "^#" *.mod.EDTA.TEanno.gff3 > no_hash.mod.EDTA.TEanno.gff3
        awk -F'\t' 'NR==1 || (\$5 - \$4) >= $threshold' no_hash.mod.EDTA.TEanno.gff3 > filtered.mod.EDTA.TEanno.gff3
        bedtools maskfasta -fi $genome -bed filtered.mod.EDTA.TEanno.gff3 -fo custom.masked.fna
        """
}



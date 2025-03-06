process GFFREAD {
    label 'gffread'
    label 'campus'

    time 2.h
    cpus 4
    memory 4.GB

    input:
        path(gtf)
        path genome

    output:
        path("*.gff"), emit: st_gff_ch
        path("*.spliced.fa"), emit: exons_ch

    script:
        """
        gffread --gtf $gtf -o ${gtf.baseName}.gff
        gffread -w ${gtf.baseName}.spliced.fa -g $genome $gtf
        """
}


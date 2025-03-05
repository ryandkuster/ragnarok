process GFFREAD {
    label 'gffread'
    label 'campus'

    time 2.h
    cpus 4
    memory 4.GB

    input:
        tuple val(sample), path(gtf)
        path genome

    output:
        tuple val(sample), path("*.gff"), emit: gff_ch
        tuple val(sample), path("*.spliced.fa"), emit: exons_ch

    script:
        """
        gffread --gtf $gtf -o ${gtf.baseName}.gff
        gffread -w ${gtf.baseName}.spliced.fa -g $genome $gtf
        """
}


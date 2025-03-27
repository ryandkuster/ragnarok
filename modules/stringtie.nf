process STRINGTIE {
    label 'stringtie'
    label 'campus'

    time 2.h
    cpus 24
    memory 80.GB

    input:
        path(bam)

    output:
        path("*.gtf"), emit: st_ch

    script:
        """
        stringtie \
            -p ${task.cpus} \
            $bam \
            -o 10kIntron_stringtie.gtf
        """
}

process STRINGTIE_MIX {
    label 'stringtie'
    label 'campus'

    time 2.h
    cpus 24
    memory 80.GB

    input:
        path(bam_ill)
        path(bam_iso)

    output:
        path("*.gtf"), emit: st_ch

    script:
        """
        stringtie \
            --mix \
            -p ${task.cpus} \
            $bam_ill \
            $bam_iso \
            -o 10kIntron_stringtie.gtf
        """
}

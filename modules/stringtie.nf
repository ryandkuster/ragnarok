process STRINGTIE {
    label 'stringtie'
    label 'campus'

    time 2.h
    cpus 24
    memory 80.GB

    input:
        tuple val(sample), path(bam)

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


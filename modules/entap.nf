process ENTAP_INI {
    label 'entap'
    label 'short'

    time 3.h
    cpus 20
    memory 20.GB

    script:
        """
        EnTAP \
            --run \
            --run-ini path/to/entap_run.params \
            --entap-ini path/to/entap_config.ini \
            -t ${task.cpus}
        """
}

process ENTAP_RUN {
    label 'entap'
    label 'campus'

    time 12.h
    cpus 20
    memory 50.GB

    input:
        path("*")

    script:
        """
        mikado.loci_out.proteins.fa \
        """
}
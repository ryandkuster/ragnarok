process ENTAP_INI {
    label 'entap'
    label 'campus'

    time 12.h
    cpus 20
    memory 60.GB

    input:
        path(entap_conf)
        path(entap_run)

    output:
        path("entap_config.ini"), emit: conf_ch
        path("entap_run.params"), emit: run_ch
        path("entap_outfiles"), emit: db_ch

    script:
        """
        cp $entap_conf entap_config.ini
        cp $entap_run entap_run.params

        EnTAP \
            --config \
            --run-ini ./entap_run.params \
            --entap-ini ./entap_config.ini \
            -t ${task.cpus}
        """
}

process ENTAP_RUN {
    label 'entap'
    label 'campus'

    publishDir(path: "${publish_dir}/entap", mode: "copy")

    time 24.h
    cpus 20
    memory { 150.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 1

    stageInMode 'copy'

    input:
        path(entap_conf)
        path(entap_run)
        path(entap_db)
        path(proteins)

    output:
        path("entap_outfiles/final_results"), emit: entap_ch
        path("entap_outfiles/final_results/annotated_without_contam.tsv"), emit: annot_ch

    script:
        """
        grep -v '^database' entap_run.params > tmp.txt
        mv tmp.txt entap_run.params
        echo "input=\$(realpath mikado.loci_out.proteins.fa)" >> entap_run.params

        mapfile -t dmnd_files < <(find ./entap_outfiles/bin/ -type f -name "*.dmnd" -exec realpath {} \\;)
        dmnd=\$(IFS=,; echo "\${dmnd_files[*]}")
        echo "database=\${dmnd}" >> entap_run.params

        EnTAP \
            --run \
            --run-ini entap_run.params \
            --entap-ini entap_config.ini \
            -t ${task.cpus}
        """
}
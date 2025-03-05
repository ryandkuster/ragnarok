process MULTIQC {
    label 'multiqc'
    label 'campus'

    time 6.h
    cpus 6
    memory 2.GB

    publishDir(path: "${publish_dir}/qc/${outdir_name}", mode: "copy")

    input:
        path('*')
        val outdir_name

    output:
        path("*html")

    script:
        """
        multiqc .
        """
}

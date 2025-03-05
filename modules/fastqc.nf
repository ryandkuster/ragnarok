process FASTQC {
    label 'fastqc'
    label 'campus'

    time 12.h
    cpus 6
    memory 2.GB

    publishDir(path: "${publish_dir}/qc/${outdir_name}", mode: "copy")

    input:
        tuple val(sample), path(r1), path(r2)
        val outdir_name

    output:
        path("fastqc_${r1.baseName}_logs"), emit: fastq_ch

    script:
        """
        mkdir fastqc_${r1.baseName}_logs
        fastqc -o fastqc_${r1.baseName}_logs -q ${r1} ${r2}
        """
}

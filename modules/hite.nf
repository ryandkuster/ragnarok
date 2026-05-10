process HITE {
    label 'hite'
    label 'long'

    time 48.h
    cpus 48
    memory { 80.GB * task.attempt }

    publishDir(path: "${params.publish_dir}/publish/masking", mode: "copy")

    errorStrategy 'retry'
    maxRetries 1

    input:
        path(genome)

    output:
        path("hite_output/HiTE.gff"), emit: hite_ch

    script:
        """
        mkdir hite_output

        python /HiTE/main.py \
            --genome $genome \
            --thread $task.cpus \
            --out_dir hite_output \
            --annotate 1 \
            --recover 1
        """
}

process FASTP_ADAPTERS {
    label "fastp"
    label "campus"

    time 12.h
    cpus 6
    memory 20.GB

    publishDir(path: "${publish_dir}/cutadapt", mode: "symlink")

    input:
        tuple val(sample), path(r1), path(r2)
        val minimum_length

    output:
        tuple val(sample), path("cut_${r1}"), path("cut_${r2}"), emit : reads

    script:
        forward = "cut_${r1}"
        reverse = "cut_${r2}"
        """
        fastp \
          --thread $task.cpus \
          -i $r1 \
          -I $r2 \
          -o $forward \
          -O $reverse \
          --length_required ${minimum_length}
        """
}

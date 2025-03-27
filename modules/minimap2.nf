process MINIMAP2 {
    label 'minimap2'
    label 'campus'

    time 10.h
    cpus 20
    memory 80.GB

    input:
        path genome
        path long_ch

    output:
        path("isoseq.sam"), emit: mp_ch

    script:
        """
        minimap2 -t "${task.cpus}" -ax splice:hq -uf "${genome}" $long_ch > isoseq.sam
        """
}


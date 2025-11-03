process MINIMAP2_PB {
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

process MINIMAP2_ONT {
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
        minimap2 -t "${task.cpus}" -ax splice -uf -k14 "${genome}" $long_ch > isoseq.sam
        """
}

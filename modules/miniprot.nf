process MINIPROT {
    label 'miniprot'
    label 'campus'

    time 10.h
    cpus 20
    memory 40.GB

    input:
        path(gff3)
        path genome
        path protein

    output:
        path("aa_miniprot.gff"), emit: mp_ch

    script:
        """
        miniprot -t${task.cpus} -d ${gff3.baseName}.mpi $genome
        miniprot -Iut${task.cpus} --gff ${gff3.baseName}.mpi $protein > aa_miniprot.gff
        """
}


process FPNLRS_SETUP {
    label 'short'

    time 3.h
    cpus 24
    memory 10.GB

    input:
        path(genome)

    output:
        path("FindPlantNLRs"), emit: fplnr_db_ch

    script:
        """
        wget https://doi.org/10.1371/journal.pbio.3001124.s013

        git clone --branch docker_version https://github.com/ZhenyanLuo/FindPlantNLRs.git
        mkdir -p FindPlantNLRs/genome
        mkdir -p FindPlantNLRs/tmp

        mv journal.pbio.3001124.s013 FindPlantNLRs/ref_db/ref.fasta
        cp $genome FindPlantNLRs/genome/test.fa
        """
}

process FINDPLANTNLRS {
    label 'findplantnlrs'
    label 'short'

    containerOptions "--bind ~ --bind $fplnr_db:/home/FindPlantNLRs/ --bind $ipscan:/home/interproscan"

    time 3.h
    cpus 40
    memory 100.GB

    input:
        path(fplnr_db)
        path(ipscan)

    output:
        path("FindPlantNLRs"), emit: fplnr_ch

    script:
        """
        /opt/conda/bin/activate FindPlantNLRs
        cd /home/FindPlantNLRs/
        snakemake \
            -s FindPlantNLRs \
            -c ${task.cpus} \
            --wait-for-files
        """
}

process ANNOTATENLRS {
    label 'findplantnlrs'
    label 'short'

    containerOptions "--bind ~ --bind $fplnr:/home/FindPlantNLRs/ --bind $ipscan:/home/interproscan --bind $genemark:/root/gmes_linux_64"

    publishDir(path: "${publish_dir}/find_plant_nlrs", mode: "copy")

    time 3.h
    cpus 40
    memory 100.GB

    input:
        path(fplnr)
        path(ipscan)
        path(genemark)

    output:
        path("FindPlantNLRs"), emit: fplnr_ch

    script:
        """
        cd /home/FindPlantNLRs
        /opt/conda/bin/activate Annotate_NLR
        export AUGUSTUS_SCRIPTS_PATH=/opt/conda/pkgs/augustus-3.5.0-pl5321heb9362c_5/bin/
        cp -r /opt/conda/pkgs/augustus-3.5.0-pl5321heb9362c_5/config/ .
        export AUGUSTUS_CONFIG_PATH=/home/FindPlantNLRs/config
        snakemake -s Annotate_NLR -c 12 --wait-for-files
        """
}

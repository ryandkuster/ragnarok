process FPNLRS_SETUP {
    label 'short'

    time 3.h
    cpus 24
    memory 10.GB

    input:
        path(genome)

    output:
        path("FindPlantNLRs"), emit: fpnlr_db_ch

    script:
        """
        wget https://doi.org/10.1371/journal.pbio.3001124.s013

        git clone --branch docker_version https://github.com/ZhenyanLuo/FindPlantNLRs.git
        mkdir -p FindPlantNLRs/genome
        mkdir -p FindPlantNLRs/tmp

        mv journal.pbio.3001124.s013 FindPlantNLRs/ref_db/ref.fasta
        cp $genome FindPlantNLRs/genome/fpnlr.fa
        """
}

process FINDPLANTNLRS {
    label 'findplantnlrs'
    label 'short'

    beforeScript 'mkdir FPNLR_workaround'
    // containerOptions "--bind ~ --bind $fpnlr_db:/home/FindPlantNLRs/ --bind $ipscan:/home/interproscan"
    containerOptions "--bind ~ --bind FPNLR_workaround:/home/FindPlantNLRs/ --bind $ipscan:/home/interproscan"

    time 3.h
    cpus 40
    memory 100.GB

    input:
        path(fpnlr_db)
        path(ipscan)

    output:
        // path("FindPlantNLRs"), emit: fpnlr_ch
        path("FPNLR_workaround"), emit: fpnlr_ch

    script:
        """
        cp -r ${fpnlr_db}/* FPNLR_workaround
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
    label 'campus'

    beforeScript 'mkdir FindPlantNLRs'

    // containerOptions "--bind ~ --bind $fpnlr:/home/FindPlantNLRs/ --bind $ipscan:/home/interproscan --bind $genemark:/root/gmes_linux_64"
    containerOptions "--bind ~ --bind FindPlantNLRs:/home/FindPlantNLRs/ --bind $ipscan:/home/interproscan --bind $genemark:/root/gmes_linux_64"

    publishDir(path: "${publish_dir}/find_plant_nlrs", mode: "copy")

    time 12.h
    cpus 40
    memory 100.GB

    input:
        path(fpnlr)
        path(ipscan)
        path(genemark)

    output:
        path("FindPlantNLRs"), emit: fpnlr_ch
        path("FindPlantNLRs/result/fpnlr_braker_aa_NLR.gff3"), emit: nlr_gff_ch

    script:
        """
        cp -r ${fpnlr}/* FindPlantNLRs
        cd /home/FindPlantNLRs
        /opt/conda/bin/activate Annotate_NLR
        export AUGUSTUS_SCRIPTS_PATH=/opt/conda/pkgs/augustus-3.5.0-pl5321heb9362c_5/bin/
        cp -r /opt/conda/pkgs/augustus-3.5.0-pl5321heb9362c_5/config/ .
        export AUGUSTUS_CONFIG_PATH=/home/FindPlantNLRs/config
        snakemake -s Annotate_NLR -c 12 --wait-for-files
        """
}

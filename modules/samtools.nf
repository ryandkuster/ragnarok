process SAM_SORT {
    label 'samtools'
    label 'campus'

    time 2.h
    cpus 24
    memory 40.GB

    publishDir(path: "${publish_dir}/alignments/sorted_bam", mode: "copy")

    input:
        path(bam_ls)
        val(prefix)

    output:
        path("*sorted_merged.bam"), emit: sort_ch

    script:
        """
        samtools merge -@ 10 -f - $bam_ls | \
        samtools sort -@ 10 -o ${prefix}_sorted_merged.bam
        """
}

process SAM_DICT {
    label 'samtools'
    label 'short'

    time 1.h
    cpus 4
    memory 2.GB

    input:
        path(genome)

    output:
        path("*.fai"), emit: fai
        path("*.dict"), emit: dict

    script:
        """
        samtools faidx $genome
        samtools dict $genome > ${genome.baseName}.dict
        """
}

process SAM_INDEX {
    label 'samtools'
    label 'short'

    time 1.h
    cpus 4
    memory 2.GB

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path(bam), path("*bai"), emit: idx_ch

    script:
        """
        samtools index \
            -@ ${task.cpus} \
            $bam
        """
}

process SAM_TABIX {
    label 'samtools'
    label 'short'

    time 1.h
    cpus 4
    memory 2.GB

    input:
        path(vcf)

    output:
        path("*.vcf.gz"), emit: vcf_ch
        path("*.tbi"), emit: tab_ch

    script:
        """
        bgzip -c $vcf > ${vcf.baseName}.vcf.gz
        tabix -p vcf ${vcf.baseName}.vcf.gz
        """
}

process SAM_STATS {
    label 'samtools'
    label 'campus'

    time 6.h
    cpus 6
    memory 2.GB

    publishDir(path: "${publish_dir}/qc/${outdir_name}", mode: "copy")

    input:
        tuple val(sample), path(bam)
        val outdir_name

    output:
        path("*stats.txt"), emit: bam_ch

    script:
        """
        samtools stats $bam > ${bam.baseName}_stats.txt
        """
}

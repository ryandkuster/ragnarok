process STAR_INDEX_NA {
    label 'star'
    label 'campus'

    time 4.h
    cpus 10
    memory 50.GB

    input:
        path(genome)

    output:
        path("Star_index"), emit: star_idx

    script:
        """
        mkdir Star_index
        mv $genome ./Star_index

        STAR \
          --runMode genomeGenerate \
          --genomeDir ./Star_index \
          --genomeSAindexNbases 13 \
          --genomeFastaFiles ./Star_index/${genome} \
          --runThreadN ${task.cpus}
        """
}

process STAR_MAP {
    label 'star'
    label 'campus'

    time 24.h
    cpus 40
    memory 50.GB

    publishDir(path: "${publish_dir}/alignments/star", mode: "symlink")

    input:
        tuple val(sample), path(r1), path(r2)
        path(star_index)
        val(max_intron)

    output:
        tuple val(sample), path("*.bam"), emit: bam_ch
        path("*out")

    script:
        """
        STAR \
          --genomeDir ${star_index} \
          --readFilesIn ${r1} ${r2} \
          --outFileNamePrefix ${sample} \
          --outSAMstrandField intronMotif \
          --runThreadN ${task.cpus} \
          --readFilesCommand zcat \
          --outSAMtype BAM Unsorted \
          --alignIntronMax $max_intron \
          >& star_${sample}.out 
        """
}

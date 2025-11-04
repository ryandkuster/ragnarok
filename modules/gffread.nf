process GFFREAD {
    label 'gffread'
    label 'campus'

    time 2.h
    cpus 4
    memory 4.GB

    input:
        path(gtf)
        path genome

    output:
        path("*.gff"), emit: st_gff_ch
        path("*.spliced.fa"), emit: exons_ch

    script:
        """
        gffread --gtf $gtf -o ${gtf.baseName}.gff
        gffread -w ${gtf.baseName}.spliced.fa -g $genome $gtf
        """
}

process GFFREAD_MIKADO {
    label 'gffread'
    label 'short'

    time 3.h
    cpus 20
    memory 20.GB

    input:
        tuple path(mk_loci), path(mk_subloci)
        path(genome)

    output:
        tuple path("mikado.loci_out.transcripts.fa"), path("mikado.loci_out.proteins.fa"), emit: final_ch

    script:
        """
        gffread \
            -w mikado.loci_out.transcripts.fa \
            -g $genome \
            $mk_loci

        gffread \
            -y mikado.loci_out.proteins.fa \
            -g $genome \
            $mk_loci
        """
}

process GFFREAD_ENTAP {
    label 'gffread'
    label 'short'

    publishDir(path: "${params.publish_dir}/publish/RAGNAROK", mode: "copy", pattern: "*.entap_filtered.*")

    time 3.h
    cpus 20
    memory 20.GB

    input:
        path(entap_gff)
        path(genome)
        val(final_prefix)

    output:
        tuple path("*.transcripts.fna"), path("*.proteins.faa"), emit: final_ch

    script:
        """
        gffread \
            -w ${final_prefix}.entap_filtered.transcripts.fna \
            -g $genome \
            $entap_gff

        gffread \
            -y ${final_prefix}.entap_filtered.proteins.faa \
            -g $genome \
            $entap_gff
        """
}
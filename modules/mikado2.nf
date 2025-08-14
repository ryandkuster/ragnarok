process MIKADO_CONF {
    label 'mikado2'
    label 'short'

    time 3.h
    cpus 20
    memory 4.GB

    publishDir(path: "${publish_dir}/mikado/mikado_in", mode: "copy")

    input:
        path(all_gffs)
        path(design)
        path(genome)
        path(scoring)
        path(homology)

    output:
        path("configuration.yaml"), emit: yaml_ch
        path("*.fai"), emit: fai_ch
        tuple path("mikado_prepared.gtf"), path("mikado_prepared.fasta"), emit: mk_ch
        path("*.{gtf,gff,gff3}"), emit: all_ch

    script:
        """
        for ext in gff gtf gff3; do
            find . -maxdepth 1 -name "*.\${ext}" -print | while read -r file; do
                file=\$(basename \$file)
                cp \$file pre_mikado_\${file}
            done
        done

        mikado configure \
            --list mikado.tsv \
            --reference $genome \
            --mode permissive \
            --scoring $scoring \
            -bt $homology \
            configuration.yaml
        mikado prepare --json-conf configuration.yaml
        """
}

process THE_GRANDMASTER {
    label 'mikado2'
    label 'short'

    time 3.h
    cpus 24
    memory 20.GB

    publishDir(path: "${publish_dir}/mikado/mikado_out", mode: "copy")

    input:
        path(yaml)
        tuple path(mk_gtf), path(mk_fasta)
        path(dmnd)
        path(orfs)
        path(homology)
        path(genome)

    output:
        tuple path("mikado.loci_out.gff3"), path("mikado.subloci.gff3"), emit: gm_ch

    script:
        """
        mikado serialise \
            --json-conf $yaml \
            --transcripts $mk_fasta \
            --xml $dmnd \
            --orfs $orfs \
            --blast_targets $homology

        mikado pick \
            -p ${task.cpus} \
            --fasta $genome \
            --configuration $yaml \
            --subloci-out mikado.subloci.gff3 \
            --loci-out mikado.loci_out.gff3
        """
}

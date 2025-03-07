process DIAMOND {
    label 'diamond'
    label 'short'

    time 3.h
    cpus 24
    memory 20.GB

    input:
        path(homology)
        tuple path(mk_gtf), path(mk_fasta)

    output:
        path("mikado_prepared.blast.tsv"), emit: dmnd_ch

    script:
        """
        diamond makedb \
            --in $homology \
            -d uniprot_sprot_diamond

        diamond blastx \
            -k 5 \
            -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop \
            -p ${task.cpus} \
            --db uniprot_sprot_diamond \
            -q $mk_fasta \
            -o mikado_prepared.blast.tsv
        """
}


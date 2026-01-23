process PSAURON {
    label 'psauron'
    label 'gpu'

    time 1.h
    cpus 2
    memory 40.GB

    publishDir(path: "${params.publish_dir}/publish/RAGNAROK", mode: "copy")

    input:
        tuple path(entap_transcripts), path(entap_proteins)
        val(final_prefix)

    output:
        path("${final_prefix}.entap_filtered.proteins.psauron.txt"), emit: ch_score

    script:
        """
        psauron \
            -p \
            -i $entap_proteins \
            -o ${final_prefix}.entap_filtered.proteins.psauron.txt
        """
}

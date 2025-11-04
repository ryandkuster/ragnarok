process AGAT_SUBSET {
    label 'agat'
    label 'short'

    publishDir(path: "${params.publish_dir}/publish/RAGNAROK", mode: "copy", pattern: "*.entap*.gff3")

    time 3.h
    cpus 20
    memory 20.GB

    input:
        tuple path(mk_loci), path(mk_subloci)
        path(entap_annot)
        val(final_prefix)

    output:
        path("*.entap_filtered.gff3"), emit: subset_pass_ch
        path("*.entap_no_annotation.gff3"), emit: subset_fail_ch

    script:
        """
        grep "gene" mikado.loci_out.gff3 | awk '{print \$9}' | cut -d "=" -f 2 | cut -d ";" -f 1 | sed -E 's/\\.[0-9]+\$//' | sort | uniq > all_ids.txt
        tail -n +2 annotated_without_contam.tsv | awk '{print \$1}' | sed -E 's/\\.[0-9]+\$//' | sort | uniq > pass_ids.txt
        comm -23 all_ids.txt pass_ids.txt > fail_ids.txt

        agat_sp_filter_feature_from_keep_list.pl \
            --gff mikado.loci_out.gff3 \
            --keep_list pass_ids.txt \
            --output ${final_prefix}.entap_filtered.gff3

        agat_sp_filter_feature_from_keep_list.pl \
            --gff mikado.loci_out.gff3 \
            --keep_list fail_ids.txt \
            --output ${final_prefix}.entap_no_annotation.gff3
        """
}
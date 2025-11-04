process PARSE_INPUT {
    label 'pandas'
    label 'short'

    time 5.m
    cpus 2
    memory 10.GB

    publishDir "${params.publish_dir}/publish/design", mode: 'copy'

    input:
        path(design)
        val(skip_st)
        val(skip_hx)
        val(nlrs)
        val(lo_genome)
        val(lo_gff)

    output:
        path 'mikado.tsv', emit: design_ch
        path 'gff_paths.csv', emit: path_ch

    script:
        """
        parse_input.py \
            $design \
            $skip_st \
            $skip_hx \
            $nlrs \
            $lo_genome \
            $lo_gff
        """
}

process PROT_FIX {
    label 'pandas'
    label 'short'

    time 12.m
    cpus 2
    memory 10.GB

    stageInMode 'copy'

    input:
        path("*")

    output:
        path 'mikado.loci_out.proteins.fa', emit: prot_ch

    script:
        """
        mv mikado.loci_out.proteins.fa original_mikado.loci_out.proteins.fa
        prot_check.py \
            original_mikado.loci_out.proteins.fa \
            mikado.loci_out.proteins.fa
        """
}

process SCAF2NUM {
    label 'pandas'
    label 'short'

    time 12.m
    cpus 2
    memory 10.GB

    input:
        path(genome)
        path(cds)

    output:
        path 'scaf2num_genome.fna', emit: s2n_ref_ch
        path 'scaf2num_cds.fna', emit: s2n_cds_ch
        path 's2n_ref.json', emit: s2n_json_ch

    script:
        """
        scaffold_to_number.py \
            $genome \
            scaf2num_genome.fna \
            s2n_ref.json

        scaffold_to_number.py \
            $cds \
            scaf2num_cds.fna \
            s2n_cds.json
        """
}

process NUM2SCAF {
    label 'pandas'
    label 'short'

    time 12.m
    cpus 2
    memory 10.GB

    input:
        path(s2n_genome)
        path(s2n_json)

    output:
        path 'num2scaf_genome.fna', emit: n2s_ref_ch

    script:
        """
        number_to_scaffold.py \
            $s2n_genome \
            num2scaf_genome.fna \
            s2n_ref.json
        """
}

process NOZIP_REF {
    label 'pandas'
    label 'short'

    time 12.m
    cpus 2
    memory 10.GB

    input:
        path(genome)

    output:
        path 'genome_unzip.fna', emit: nozip_ref_ch

    script:
        """
        if [[ \$(head -c 2 "$genome" | od -An -tx1) =~ "1f 8b" ]]; then
            cp $genome genome_unzip.fna.gz
            gunzip genome_unzip.fna.gz
        else
            cp $genome genome_unzip.fna
        fi
        """
}

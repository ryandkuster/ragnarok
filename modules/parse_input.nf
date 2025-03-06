process PARSE_INPUT {
    label 'pandas'
    label 'short'

    time 5.m
    cpus 2
    memory 1.GB

    publishDir "${publish_dir}/design", mode: 'copy'

    input:
        path(design)
        val(skip_st)
        val(skip_hx)

    output:
        path 'mikado.txt', emit: design_ch

    script:
        """
        parse_input.py \
            $design \
            $skip_st \
            $skip_hx
        """
}

process COMPLEASM_DB {
    label 'compleasm'
    label 'short'
  
    time 3.h
    cpus 20
    memory 4.GB

    input:
        val(busco_db)

    output:
        path("mb_downloads"), emit: db_ch
  
    script:
        """
        compleasm download $busco_db
        """
}

process COMPLEASM {
    label 'compleasm'
    label 'short'
  
    time 3.h
    cpus 20
    memory 40.GB

    publishDir(path: "${publish_dir}/compleasm", mode: "copy")
  
    input:
        path("*")
        path("*")
        val(busco_db)
  
    output:
        path("compleasm/*"), emit: ch_score
  
    script:
        """
        compleasm protein \
            -p mikado.loci_out.proteins.fa \
            -l $busco_db \
            -t ${task.cpus} \
            -o compleasm 
        """
}
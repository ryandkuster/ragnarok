process COMPLEASM_DB {
    label 'compleasm'
    label 'short'
  
    time 3.h
    cpus 20
    memory 4.GB
  
    output:
        path("mb_downloads"), emit: db_ch
  
    script:
        """
        compleasm download embryophyta_odb10
        """
}

process COMPLEASM {
    label 'compleasm'
    label 'short'
  
    time 3.h
    cpus 20
    memory 4.GB

    publishDir(path: "${publish_dir}/compleasm", mode: "copy")
  
    input:
        path("*")
        path("*")
  
    output:
        path("compleasm/*"), emit: ch_score
  
    script:
        """
        compleasm protein \
            -p mikado.loci_out.proteins.fa \
            -l embryophyta_odb10 \
            -t ${task.cpus} \
            -o compleasm 
        """
}
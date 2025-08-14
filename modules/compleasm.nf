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

    publishDir(path: "${publish_dir}/RAGNAROK", mode: "copy")
  
    input:
        tuple path(entap_transcripts), path(entap_proteins)
        path("*")
        val(busco_db)
        val(final_prefix)
  
    output:
        path("${final_prefix}.entap_filtered.compleasm_${busco_db}.txt"), emit: ch_score
  
    script:
        """
        compleasm protein \
            -p $entap_proteins \
            -l $busco_db \
            -t ${task.cpus} \
            -o compleasm
        
        cp compleasm/summary.txt ${final_prefix}.entap_filtered.compleasm_${busco_db}.txt
        """
}
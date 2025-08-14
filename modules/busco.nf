process BUSCO {
    label 'busco'
    label 'short'
  
    time 3.h
    cpus 20
    memory 40.GB

    publishDir(path: "${publish_dir}/RAGNAROK", mode: "copy")
  
    input:
        tuple path(entap_transcripts), path(entap_proteins)
        val(busco_db)
        val(final_prefix)
  
    output:
        path("${final_prefix}.entap_filtered.busco_${busco_db}.txt"), emit: ch_score
  
    script:
        """
        busco \
            -m transcriptome \
            -i $entap_transcripts \
            -c ${task.cpus} \
            -o busco_transcriptome \
            -l $busco_db
        
        cp busco_transcriptome/short_summary*.txt ${final_prefix}.entap_filtered.busco_${busco_db}.txt
        """
}

process BUSCO {
    label 'busco'
    label 'short'
  
    time 3.h
    cpus 20
    memory 40.GB

    publishDir(path: "${publish_dir}/busco", mode: "copy")
  
    input:
        path("*")
  
    output:
        path("busco_transcriptome/*"), emit: ch_score
  
    script:
        """
        busco \
            -m transcriptome \
            -i ./mikado.loci_out.transcripts.fa \
            -c ${task.cpus} \
            -o busco_transcriptome \
            -l embryophyta_odb10
        """
}

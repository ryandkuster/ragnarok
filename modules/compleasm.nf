process COMPLEASM {
    label 'compleasm'
    label 'short'
  
    time 3.h
    cpus 20
    memory 4.GB

    publishDir(path: "${publish_dir}/compleasm", mode: "copy")
  
    input:
        path("*")
  
    output:
        path("compleasm/*"), emit: ch_score
  
    script:
        """
        compleasm.py protein \
            -p mikado.loci_out.proteins.fa \
            -l embryophyta_odb10 
            -t ${task.cpus} \
            -o compleasm 
        """
}

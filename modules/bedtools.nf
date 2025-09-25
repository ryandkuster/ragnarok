process EDTA_THRESHOLD {
    label 'bedtools'
    label 'short'

    time 3.h
    cpus 24
    memory 10.GB

    publishDir(path: "${publish_dir}/edta_masking", mode: "copy")

    input:
        path("*")
        path(genome)
        val(threshold)

    output:
        path("custom.masked.fna"), emit: threshold_ch

    script:
        """
        grep -v "^#" *.mod.EDTA.TEanno.gff3 > no_hash.mod.EDTA.TEanno.gff3
        awk -F'\t' 'NR==1 || (\$5 - \$4) >= $threshold' no_hash.mod.EDTA.TEanno.gff3 > filtered.mod.EDTA.TEanno.gff3
        bedtools maskfasta -fi $genome -bed filtered.mod.EDTA.TEanno.gff3 -fo custom.masked.fna
        """
}

process HITE_THRESHOLD {
    label 'bedtools'
    label 'short'

    time 3.h
    cpus 24
    memory 10.GB

    publishDir(path: "${publish_dir}/edta_masking", mode: "copy")

    input:
        path("*")
        path(genome)
        val(threshold)

    output:
        path("custom.masked.fna"), emit: threshold_ch

    script:
        """
        grep -v "^#" HiTE.gff > no_hash.HiTE.gff3
        awk -F'\t' 'NR==1 || (\$5 - \$4) >= $threshold' no_hash.HiTE.gff3 > filtered.HiTE.gff3
        bedtools maskfasta -fi $genome -bed filtered.HiTE.gff3 -fo custom.masked.fna
        """
}

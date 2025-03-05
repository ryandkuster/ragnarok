process CONCAT {
    label 'bcftools'
    label 'campus'

    time 12.h
    cpus 2

    input:
        path vcfs

    output:
        path("combined_intervals.vcf"), emit: merge_ch

    script:
        """
        inputs=""
        for i in \$( ls -v *combined.vcf ) ; do
            inputs+="\$i "
        done

        bcftools concat -o combined_intervals.vcf \${inputs}
        """
}

process CONCAT_GZ {
    label 'bcftools'
    label 'campus'

    time 12.h
    cpus 2

    input:
        path vcfs

    output:
        path("combined_intervals.vcf.gz"), emit: merge_ch

    script:
        """
        inputs=""
        for i in \$( ls -v *.vcf.gz ) ; do
            inputs+="\$i "
        done

        bcftools concat -o combined_intervals.vcf.gz \${inputs}
        """
}

process BCF_SORT {
    label 'bcftools'
    label 'long'

    time 48.h
    cpus 20
    memory 100.GB

    input:
        path vcf

    output:
        path("${vcf}.gz"), emit: bcf_sort_ch

    script:
        """
        mkdir tmp
        bcftools sort --temp-dir ./tmp $vcf -Oz -o ${vcf}.gz
        """
}

process BCF_IDX {
    label 'bcftools'
    label 'campus'

    time 24.h
    cpus 20
    memory 100.GB

    input:
        path(vcf)

    output:
        path("*tbi"), emit: vcf_ch

    script:
        """
        bcftools index -t -f ${vcf}
        """
}

process BCF_STATS {
    label 'bcftools'
    label 'campus'

    time 24.h
    cpus 20
    memory 100.GB

    input:
        path(vcf)
        path(tbi)
        val(samples)

    output:
        path("bcftools_stats.txt"), emit: bcf_stats_ch

    script:
        inputs = ""
        for (i in samples) {
          inputs += "${i},"
        }

        """
        bcftools stats ${vcf} -s $inputs > bcftools_stats.txt
        """
}

process BCF_MPILEUP {
    label 'bcftools'
    label 'campus'

    time 24.h
    cpus 20
    memory 50.GB

    input:
        path(genome)
        path(region)
        path(bam)

    output:
        path("*vcf.gz"), emit: vcf_ch

    script:
        """
        BAMS=\$(ls -1 *bam)
        echo \$BAMS

        bcftools mpileup -Ou \
            --annotate FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/SCR \
            -f $genome \
            -R $region \$BAMS | \
            bcftools call -mv -Oz -o ${region.baseName}.vcf.gz
        """
}

process BCF_SNP {
    label 'bcftools'
    label 'campus'

    publishDir(path: "${publish_dir}/${outdir_name}/variants", mode: "copy")

    time 12.h
    cpus 2

    input:
        path vcf
        path idx
        val filt_params
        val(outdir_name)

    output:
        path("*_final*"), emit: final_ch

    script:
        """
        bcftools view \
            -v snps \
            -i '$filt_params' \
            $vcf \
            -O z \
            -o snp_final.vcf.gz
        """
}

process BCF_IND {
    label 'bcftools'
    label 'campus'

    publishDir(path: "${publish_dir}/${outdir_name}/variants", mode: "copy")

    time 12.h
    cpus 2

    input:
        path vcf
        path idx
        val filt_params
        val(outdir_name)

    output:
        path("*_final*"), emit: final_ch

    script:
        """
        bcftools view \
            -v indels \
            -i '$filt_params' \
            $vcf \
            -O z \
            -o indel_final.vcf.gz
        """
}

process BCF_GATK_SNP {
    label 'bcftools'
    label 'campus'

    publishDir(path: "${publish_dir}/${outdir_name}/variants", mode: "copy")

    time 12.h
    cpus 2

    input:
        path vcf
        path idx
        val filt_params
        val(outdir_name)

    output:
        path("*_final*"), emit: final_ch

    script:
        """
        bcftools view \
            -i '$filt_params' \
            -f PASS \
            $vcf \
            -O z \
            -o snp_final.vcf.gz
        """
}

process BCF_GATK_IND {
    label 'bcftools'
    label 'campus'

    publishDir(path: "${publish_dir}/${outdir_name}/variants", mode: "copy")

    time 12.h
    cpus 2

    input:
        path vcf
        path idx
        val filt_params
        val(outdir_name)

    output:
        path("*_final*"), emit: final_ch

    script:
        """
        bcftools view \
            -i '$filt_params' \
            -f PASS \
            $vcf \
            -O z \
            -o indel_final.vcf.gz
        """
}

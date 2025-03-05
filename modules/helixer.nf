process HELIXER {
    label 'helixer'
    label 'gpu'

    time 4.h
    cpus 24
    memory 20.GB

    input:
        path(genome)
        val(lineage)

    output:
        path("*.gff3"), emit: gff_ch

    script:
        """
        wget https://zenodo.org/records/10836346/files/land_plant_v0.3_a_0080.h5?download=1
        Helixer.py \
            --model-filepath land_plant_v0.3_a_0080.h5 \
            --fasta-path $genome \
            --species helixer_species \
            --gff-output-path helixer.gff3
        """
}


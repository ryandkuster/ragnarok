#!/bin/bash

################################################################################
#                                   RAGNAROK!                                  #
#                         RApid Genome anNotAtion ROcKs!                       #
#                          Chris Gottschalk 12/6/2024                          #
#                                    v1.0                                      #
#                            Example with M. ionesis                           #
################################################################################
#                       0. Reading command-line input                          #
################################################################################

help_flag=false
# Optional
output='ragnarok_out'
t_num='1'

print_usage() {
        echo "Usage: $0 <-f> {GENOME_FASTA} [-h|-o {FILENAME PREFIX}|-t {NUM THREADS}]"
        echo
        echo "To display this screen use -h for help"
        echo
        echo "Required:"
        echo "-f        Masked genome fasta"
        echo "-m        Method used to generate RNA file(s). iso will use minimap2 for PacBio IsoSeq reads. ill will use TBH for Illumina reads"
        echo "-r        File directory of RNA fastas or list of file names or directory names, separate names with spaces"
        echo
        echo "Optional:"
        echo "-t        Number of threads to use, DEFAULT: 1"
        echo "-o        Output directory name, DEFAULT: ./ragnarok_out"
}

while getopts ':ho:t:f:' flag
do
        case "${flag}" in
        f) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -f requires an argument" && print_usage && exit 1; fi
           masked_fasta="${OPTARG}" ;;
        m) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -m requires input argument" && print_usage && exit 1; fi
           if [ "$OPTARG" != "iso" ] && [ "$OPTARG" != "ill" ]; then echo "Error: -m must be iso or ill" && print_usage && exit 1; fi
           rna_type="${OPTARG}" ;;
        h) help_flag=true ;;
        t) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -t requires an argument" && print_usage && exit 1; fi
           t_num="${OPTARG}" ;;
        o) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "_" ]; then echo "Error: -o requires an argument" && print_usage && exit 1; fi
           output="${OPTARG}" ;;
        r) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]; then echo "Error: -r requires input argument" && print_usage && exit 1; fi
           rna_seq+=("${OPTARG}")
           while [[ ${!OPTIND} && ${!OPTIND} != -* ]]; do rna_seq+=("${!OPTIND}") && ((OPTIND++)); done ;;
        \?) echo "Error: Invalid option -$OPTARG"
            print_usage
            exit 1 ;;
        :)  echo "option -$OPTARG requires an argument"
            print_usage
            exit 1 ;;
        *) print_usage
           exit 1;;
        esac
done

if [ "$help_flag" == true ]; then print_usage && exit 1; fi
if [ -z "$rna_files" ]; then echo "Error: No RNA files input or found!" && print_usage && exit 1; fi

mkdir -p "$output"
cd "$output"

################################################################################
#                          Step 1: Assembling Annotations                      #
################################################################################

################################################################################
#                             Helixer ab initio                                #
################################################################################

# Run the helixer.sh script, following it's instructions. Since this is the only
# step requiring gpus, it gets it's own job. 

################################################################################
#                   StringTie2 Transcriptome Assembly                          #
################################################################################

spec_name=$(echo "${masked_fasta}" | sed 's/.fasta.mod.MAKER.masked//')

echo "Activating minimap2 conda"
conda activate minimap2

#long-read isoseq option
if [ "$read_type" == "iso" ]
  then 
    minimap2 -t "${t_num}" -ax splice:hq -uf "${masked_fasta}" "${rna_seq}" > isoseq.sam
    samtools view -bS --threads "${t_num}" isoseq.sam -o isoseq.bam
    rm isoseq.sam #saving space on disk drive

    #need to remove "None" typing mappings that occur with minimap2
    samtools view -@ "${t_num}" -b -F 4 isoseq.bam > isoseq_filtered.bam
    samtools sort -@ "${t_num}" isoseq_filtered.bam -T isoseq_filtered_sort -o isoseq_filtered_sort.bam
    samtools index -@ "${t_num}" isoseq_filtered_sort.bam
    stringtie -p "${t_num}"  -L isoseq_filtered_sort.bam -o stringtie_LR.gtf
    gffread stringtie_LR.gtf stringtie_LR.gff
#short read RNAseq option
elif [ "$read_type" == "ill" ]
  then
    bamfiles=$(find "${rna_files[@]}" -type f)
    STAR --runMode genomeGenerate --runThreadN "${t_num}" --genomeDir ./fasta.mod.MAKER.masked.index --genomeFastaFiles "${masked_fasta}"
    STAR --genomeDir ./fasta.mod.MAKER.masked.index --outFileNamePrefix "${spec_name}" --readFilesCommand zcat --runThreadN "${t_num}" --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --alignIntronMax 10000 \
    --readFilesIn "${bamfiles}"

    samtools sort -o "${spec_name}"Aligned.out.sort.bam -@ "${t_num}" "${spec_name}"Aligned.out.bam
    samtools index -@ "${t_num}" "${spec_name}"Aligned.out.sort.bam
    stringtie -p "${t_num}" "${spec_name}"Aligned.out.sort.bam -o "${spec_name}"_10kIntron_stringtie.gtf
    gffread "${spec_name}"_10kIntron_stringtie.gtf "${spec_name}"_10kIntron_stringtie.gff
    gffread -w "${spec_name}"_10kIntron_stringtie.fa -g "${masked_fasta}" "${spec_name}"_10kIntron_stringtie.gtf
fi

echo "Exiting minimap2 conda"
conda deactivate

#Generate a EST based longest ORF annotation using the StringTie2 transcripts
# this might be in util, but try anyway
apptainer exec -e /project/pomics/programs/transdecoder.v5.7.1.simg gtf_genome_to_cdna_fasta.pl stringtie_LR.gtf "${masked_fasta}" > transcripts.fasta
apptainer exec -e /project/pomics/programs/transdecoder.v5.7.1.simg TransDecoder.LongOrfs -t transcripts.fasta
apptainer exec -e /project/pomics/programs/transdecoder.v5.7.1.simg TransDecoder.Predict -t transcripts.fasta
apptainer exec -e /project/pomics/programs/transdecoder.v5.7.1.simg gtf_to_alignment_gff3.pl stringtie_LR.gtf > transcripts.gff3
apptainer exec -e /project/pomics/programs/transdecoder.v5.7.1.simg cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     transcripts.gff3 \
     transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3

################################################################################
#                   Miniprot protein2genome mapping annotation                 #
################################################################################
#download or get the M.fusca hap1 protein file as fusca has the best annotation IMO

echo "Activating busco conda"
conda activate busco #has miniprot in it
miniprot -t16 -d M_ioensis_hap1_all_scaffolds.fasta.mod.MAKER.masked.mpi ./M_ioensis_hap1_all_scaffolds.fasta.mod.MAKER.masked
miniprot -Iut16 --gff M_ioensis_hap1_all_scaffolds.fasta.mod.MAKER.masked.mpi /media/chris/Drive_2/genomes/fusca/final_asm/Hap1/Mfusca_v1.0_hap1.proteins.fa > fusca_aa_miniprot.gff

echo "Exiting busco conda"
conda deactivate

################################################################################
#                       Mikado2 to select best transcripts                     #
################################################################################

conda activate mikado2

nano list.txt #need to copy and paste the text below for the list file
  """
  M_ioensis_hap1_helixer_rd1.gff3	hx	True		False	False
  stringtie_LR.gtf	st	True	1	False	True
  transcripts.fasta.transdecoder.genome.gff3	tr	False	-0.5	False	False
  fusca_aa_miniprot.gff	mp	True	1	False	False
  """

mikado configure --list list.txt --reference M_ioensis_hap1_all_scaffolds.fasta.mod.MAKER.masked --mode permissive --scoring plants.yaml -bt /home/chris/Documents/databases/uniprot_sprot.fasta configuration.yaml
mikado prepare --json-conf configuration.yaml

conda deactivate

#apply functional annotation from uniprot/swissprot database
conda activate diamond

diamond makedb --in ~/Documents/databases/uniprot_sprot.fasta -d uniprot_sprot_diamond
diamond blastx -k 5 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop -p 32 --db ~/Documents/databases/uniprot_sprot_diamond -o mikado_prepared.blast.tsv -q mikado_prepared.fasta

conda deactivate

#going to select the longsest ORFs from the Mikado prepared transcripts
conda activate transdecoder

TransDecoder.LongOrfs -t mikado_prepared.fasta

conda deactivate

conda activate mikado2

mikado serialise --json-conf configuration.yaml --transcripts mikado_prepared.fasta --xml mikado_prepared.blast.tsv --orfs mikado_prepared.fasta.transdecoder.gff3 --blast_targets ~/Documents/databases/uniprot_sprot.fasta
mikado pick -p 32 --fasta M_ioensis_hap1_all_scaffolds.fasta.mod.MAKER.masked --configuration configuration.yaml --subloci-out mikado.subloci.gff3 --loci-out mikado.loci_out.gff3

conda deactivate


################################################################################
#                       Generating final files and QC                          #
################################################################################

gffread -w mikado.loci_out.transcripts.fa -g M_ioensis_hap1_all_scaffolds.fasta.mod.MAKER.masked mikado.loci_out.gff3
gffread -y mikado.loci_out.proteins.fa -g M_ioensis_hap1_all_scaffolds.fasta.mod.MAKER.masked mikado.loci_out.gff3

conda activate busco

busco -m transcriptome -i ./mikado.loci_out.transcripts.fa -c 32 -o busco_transcriptome -l embryophyta_odb10 &
---------------------------------------------------
|Results from dataset embryophyta_odb10            |
---------------------------------------------------
|C:98.6%[S:58.4%,D:40.3%],F:0.6%,M:0.8%,n:1614     |
|1592    Complete BUSCOs (C)                       |
|942    Complete and single-copy BUSCOs (S)        |
|650    Complete and duplicated BUSCOs (D)         |
|9    Fragmented BUSCOs (F)                        |
|13    Missing BUSCOs (M)                          |
|1614    Total BUSCO groups searched               |
---------------------------------------------------

agat_sp_statistics.pl --gff mikado.loci_out.gff3 -o annotation.stats

#note Mikado will annotation ncRNA_genes. Subtract this number from the annotation.stats
#gene number to get a final protein coding number

#52,554 total genes annotated
#7,161 ncRNA_genes annotated
#=45,393 coding genes annotated

conda deactivate

#do some renaming of genes
conda activate maker

maker_map_ids --prefix drMalIon.v2.0. --justify 8 --iterate 1 mikado.loci_out.gff3 > id_map
cp mikado.loci_out.gff3 drMalIon.v2.0.genes.gff3
cp mikado.loci_out.transcripts.fa drMalIon.v2.0.transcripts.fa
cp mikado.loci_out.proteins.fa drMalIon.v2.0.proteins.fa

map_gff_ids id_map drMalIon.v2.0.genes.gff3
map_fasta_ids id_map drMalIon.v2.0.proteins.fa
map_fasta_ids id_map drMalIon.v2.0.transcripts.fa

conda deactivate

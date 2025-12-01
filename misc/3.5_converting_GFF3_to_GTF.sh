#!/bin/bash

#SBATCH --job-name=3.5_converting_GFF3_to_GTF
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A_%a.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A_%a.err

module load Anaconda3
source activate cufflinks

# Convert exon entries to pseudo-CDS
awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon"{$3="CDS"}1' /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/data/P_pacificus/Hybrid1/genes.gtf > \
						   /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/data/P_pacificus/Hybrid1/genes.gff3

# Extract CDS sequences from gff3
gffread /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/data/P_pacificus/Hybrid1/genes.gff3 -g \
        /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/data/P_pacificus/Hybrid1/sequences.fa -x \
        /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/data/P_pacificus/Hybrid1/cds.fa

# Convert gff3 to gtf
gffread /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/data/P_pacificus/Hybrid1/genes.gff3 -T -o /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/data/P_pacificus/Hybrid1/genes.gtf



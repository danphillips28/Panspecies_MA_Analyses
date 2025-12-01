#!/bin/bash
#SBATCH --job-name=3_renaming_Ncrassa_contigs.sh   
#SBATCH --output=/home/ocdm0351/DPhil/logs/3_renaming_Ncrassa_contigs_%j.log  # where %j will be replaced by the job ID
#SBATCH --error=/home/ocdm0351/DPhil/logs/3_renaming_Ncrassa_contigs_%j.err   # optional: separate error file

sed -i \
-e 's/Supercontig_12.1/I/g' \
-e 's/Supercontig_12.2/II/g' \
-e 's/Supercontig_12.3/III/g' \
-e 's/Supercontig_12.4/IV/g' \
-e 's/Supercontig_12.5/V/g' \
-e 's/Supercontig_12.6/VI/g' \
-e 's/Supercontig_12.7/VII/g' \
-e 's/Supercontig_12.8/MtDNA/g' \
/home/ocdm0351/DPhil/Groomed_VCF_DATA/Villalbadelapena_2023_MA_Ncrassa_Variants_Simple_Groomed.vcf

# Keep only canonical chromosomes
awk 'BEGIN {OFS="\t"} /^#/ {print; next} $1 ~ /^(I|II|III|IV|V|VI|VII|MtDNA)$/' \
/home/ocdm0351/DPhil/Groomed_VCF_DATA/Villalbadelapena_2023_MA_Ncrassa_Variants_Simple_Groomed.vcf \
> tmp && mv tmp /home/ocdm0351/DPhil/Groomed_VCF_DATA/Villalbadelapena_2023_MA_Ncrassa_Variants_Simple_Groomed.vcf

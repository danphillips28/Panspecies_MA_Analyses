#!/bin/bash
#SBATCH --job-name=2_renaming_Goldmann_chromosomes.sh   
#SBATCH --output=/home/ocdm0351/DPhil/logs/2_renaming_Goldmann_chromosomes_%j.log  # where %j will be replaced by the job ID
#SBATCH --error=/home/ocdm0351/DPhil/logs/2_renaming_Goldmann_chromosomes_%j.err   # optional: separate error file

sed -i -e 's/chr//g' /home/ocdm0351/DPhil/Groomed_VCF_DATA/Goldmann_2016_Trio_Hsapiens_Variants_Simple_Groomed.vcf

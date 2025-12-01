#!/bin/bash

#SBATCH --job-name=1_Removing_Bad_Samples_from_Behringer.sh
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%j.out   # Standard output (%x = job name, %j = job ID)
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%j.err    # Standard error

module load BCFtools

cd /home/ocdm0351/DPhil/Raw_Data/

bcftools view --force-samples -s ^Sample9,Sample20,Sample25,Sample33,Sample34,Sample36,Sample41,Sample43,Sample44,Sample47,Sample48,Sample57,Sample61,Sample66,Sample86,Sample95,Sample96 Behringer_2016_MA_Spombe_Variants1.vcf \
  | bcftools view -e 'AC=0' \
  > Behringer_2016_MA_Spombe_Variants1_clean.tsv

bcftools view --force-samples -s ^Sample9,Sample20,Sample25,Sample33,Sample34,Sample36,Sample41,Sample43,Sample44,Sample47,Sample48,Sample57,Sample61,Sample66,Sample86,Sample95,Sample96 Behringer_2016_MA_Spombe_Variants2.vcf \
  | bcftools view -e 'AC=0' \
  > Behringer_2016_MA_Spombe_Variants2_clean.tsv

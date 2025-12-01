#!/usr/bin/env bash
set -euo pipefail

module load Anaconda3
source activate UCSC_liftOver

# Config
OUTDIR="/home/ocdm0351/DPhil/R_Data/"
URL="https://www.pnas.org/doi/suppl/10.1073/pnas.1608162113/suppl_file/pnas.1608162113.sd01.xlsx"

cd /home/ocdm0351/DPhil/R_Data/

# Convert csv to tsv
awk 'BEGIN{FS=","; OFS="\t"} {$1="chr"$1; $1=$1; print}' \
  C_elegans_Chromatin_Domains_prelift.tsv > C_elegans_Chromatin_Domains_prelift2.tsv

# LiftOver 
echo "Running UCSC liftOver: liftOver on Celegans Chromatin Domains"
liftOver C_elegans_Chromatin_Domains_prelift2.tsv \
        /home/ocdm0351/DPhil/liftingOver/chain_files/C_elegans/ce10ToCe11.over.chain \
        C_elegans_Chromatin_Domains.tsv \
        /home/ocdm0351/DPhil/logs/unmapped_C_elegans_Chromatin_Domains.tsv

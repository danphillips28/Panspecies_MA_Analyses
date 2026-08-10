#!/bin/bash

#SBATCH --job-name=16_MutationRate_Functional_Enrichment
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%j.out
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%j.err

module load Anaconda3
source activate pandoc_env

# optional debug
which R
which Rscript

RMD_FILE="/home/ocdm0351/DPhil/scripts/MutationRate_Functional_Enrichment.Rmd"
OUTPUT_DIR="/home/ocdm0351/DPhil/R_Data/htmls"

Rscript -e "rmarkdown::render('$RMD_FILE', output_dir = '$OUTPUT_DIR')"

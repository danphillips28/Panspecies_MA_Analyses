#!/bin/bash

#SBATCH --job-name=20_Associating_Epigenetics_and_Expression_etc_with_Mutation_Rate_Submitter.sh
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%j.out
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%j.err
#SBATCH --partition=himem-gen

module load Anaconda3
source activate pandoc_env

# 🔥 THIS IS THE FIX
export R_LIBS_USER=/home/ocdm0351/R/x86_64-conda-linux-gnu-library/4.5

# Optional debug (very useful)
echo "R location:"
which R
echo "Rscript location:"
which Rscript

Rscript -e "print(.libPaths())"

RMD_FILE="/home/ocdm0351/DPhil/scripts/Associating_MutationRate_with_Epigenetics_ExpressionAgain_and_Covariates.Rmd"
OUTPUT_DIR="/home/ocdm0351/DPhil/R_Data/htmls"
Rscript -e "rmarkdown::render('$RMD_FILE', output_dir = '$OUTPUT_DIR')"





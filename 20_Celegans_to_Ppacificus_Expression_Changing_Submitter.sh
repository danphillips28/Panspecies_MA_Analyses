#!/bin/bash

#SBATCH --job-name=20_Celegans_to_Ppacificus_Expression_Changing_Submitter
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%j.out   # Standard output (%x = job name, %j = job ID)
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%j.err    # Standard error

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

RMD_FILE="/home/ocdm0351/DPhil/scripts/Celegans_to_Ppacificus_Expression_Changing.Rmd"
OUTPUT_DIR="/home/ocdm0351/DPhil/R_Data/htmls"
Rscript -e "rmarkdown::render('$RMD_FILE', output_dir = '$OUTPUT_DIR')"




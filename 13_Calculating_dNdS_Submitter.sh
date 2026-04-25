#!/bin/bash

#SBATCH --job-name=12_CalculatingDNDS_Submitter
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%j.out   # Standard output (%x = job name, %j = job ID)
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%j.err    # Standard error

module load R/4.4.1-gfbf-2023b
module load Anaconda3
source activate pandoc_env

export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4

RMD_FILE="/home/ocdm0351/DPhil/scripts/CalculatingDNDS.Rmd"
OUTPUT_DIR="/home/ocdm0351/DPhil/R_Data/htmls"
Rscript -e "rmarkdown::render('$RMD_FILE', output_dir = '$OUTPUT_DIR')"



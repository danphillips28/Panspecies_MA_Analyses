#!/bin/bash

#SBATCH --job-name=17_Downloading_P_pacificus_TimeSeries_Expression_Submitter
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%j.out   # Standard output (%x = job name, %j = job ID)
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%j.err    # Standard error

module load R/4.4.1-gfbf-2023b
module load Anaconda3
source activate pandoc_env

RMD_FILE="/home/ocdm0351/DPhil/scripts/Downloading_P_pacificus_TimeSeries_Expression.Rmd"
OUTPUT_DIR="/home/ocdm0351/DPhil/R_Data/htmls"
Rscript -e "rmarkdown::render('$RMD_FILE', output_dir = '$OUTPUT_DIR')"

cp /home/ocdm0351/DPhil/R_Data/P_pacificus_Embryonic_TimeSeries_12_Expression.tsv /home/ocdm0351/DPhil/R_Data/P_pacificus_Embryonic_Expression.tsv

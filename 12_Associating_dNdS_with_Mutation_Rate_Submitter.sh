#!/bin/bash

#SBATCH --job-name=12_Associating_dNdS_with_Mutation_Rate_Submitter
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%j.out   # Standard output (%x = job name, %j = job ID)
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%j.err    # Standard error

module load CMake/3.18.4
module load R/4.4.1-gfbf-2023b
module load Anaconda3
source activate pandoc_env

RMD_FILE="/home/ocdm0351/DPhil/scripts/Associating_dNdS_with_Mutation_Rate.Rmd"
OUTPUT_DIR="/home/ocdm0351/DPhil/R_Data/htmls"
Rscript -e "rmarkdown::render('$RMD_FILE', output_dir = '$OUTPUT_DIR')"


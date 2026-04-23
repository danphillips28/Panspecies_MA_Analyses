#!/bin/bash -l
#SBATCH --job-name=13_Associating_dNdS_with_Mutation_Rate_Submitter
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%j.out
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%j.err
#SBATCH --partition himem-gen

module purge
module load Anaconda3
source activate pandoc_env
module load R

# Force single-threaded BLAS for debugging (you already have this; keep if useful)
export OPENBLAS_CORETYPE=generic
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Using Rscript: $(which Rscript)"

# Render the R Markdown directly (safer than running a purlled script)
#Rscript --vanilla -e 'rmarkdown::render("/home/ocdm0351/DPhil/scripts/Associating_dNdS_with_Synonymous_Mutation_Rate.Rmd", output_dir = "/home/ocdm0351/DPhil/R_Data/htmls", quiet = FALSE)'
Rscript --vanilla -e 'rmarkdown::render("/home/ocdm0351/DPhil/scripts/Associating_dNdS_with_Mutation_Rate.Rmd", output_dir = "/home/ocdm0351/DPhil/R_Data/htmls", quiet = FALSE)'
#Rscript --vanilla -e 'rmarkdown::render("/home/ocdm0351/DPhil/scripts/Making_Overall_dNdS_vs_Mutation_Rate_Plot.Rmd", output_dir = "/home/ocdm0351/DPhil/R_Data/htmls", quiet = FALSE)'


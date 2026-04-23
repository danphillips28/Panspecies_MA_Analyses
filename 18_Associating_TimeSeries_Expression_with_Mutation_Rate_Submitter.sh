#!/bin/bash -l
#SBATCH --job-name=18_Associating_TimeSeries_Expression_with_Mutation_Rate_Submitter
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%j.out
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%j.err
#SBATCH --partition himem-gen

module purge
module load Anaconda3
source activate pandoc_env
module load R

echo "Using Rscript: $(which Rscript)"

# Render the R Markdown directly (safer than running a purlled script)
Rscript --vanilla -e 'rmarkdown::render("/home/ocdm0351/DPhil/scripts/Associating_MutationRate_with_TimeSeries_Expression.Rmd", output_dir = "/home/ocdm0351/DPhil/R_Data/htmls", quiet = FALSE)'


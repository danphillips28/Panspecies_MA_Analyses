#!/bin/bash -l
#SBATCH --job-name=14_Associating_dNdS_with_Mutation_Rate_Submitter
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%j.out
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%j.err
#SBATCH --partition=himem-gen
#SBATCH --mem=128G

# ----------------------------
# Environment setup
# ----------------------------
module load Miniconda3/25.7.0-2
source /apps/software/Miniconda3/25.7.0-2/etc/profile.d/conda.sh
source activate /home/ocdm0351/.conda/envs/pandoc_env

# Thread control (good practice on HPC)
export OPENBLAS_CORETYPE=generic
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# ----------------------------
# Debug info (VERY useful)
# ----------------------------
echo "Using Rscript: $(which Rscript)"
echo "Using R: $(which R)"
echo "Using pandoc: $(which pandoc)"

Rscript -e 'cat("R version:", R.version.string, "\n")'
Rscript -e 'cat("Library paths:\n"); print(.libPaths())'

# ----------------------------
# Run RMarkdown scripts
# ----------------------------

#Rscript -e 'rmarkdown::render(
#  "/home/ocdm0351/DPhil/scripts/Associating_dNdS_with_Mutation_Rate.Rmd",
#  output_dir = "/home/ocdm0351/DPhil/R_Data/htmls",
#  quiet = FALSE
#)'

#Rscript -e 'rmarkdown::render(
#  "/home/ocdm0351/DPhil/scripts/Associating_dNdS_with_Synonymous_Mutation_Rate.Rmd",
#  output_dir = "/home/ocdm0351/DPhil/R_Data/htmls",
#  quiet = FALSE
#)'

#Rscript -e 'rmarkdown::render(
#  "/home/ocdm0351/DPhil/scripts/Associating_dNdS_with_Intronic_Mutation_Rate.Rmd",
#  output_dir = "/home/ocdm0351/DPhil/R_Data/htmls",
#  quiet = FALSE
#)'

Rscript -e 'rmarkdown::render(
  "/home/ocdm0351/DPhil/scripts/Making_Overall_dNdS_vs_Mutation_Rate_Plot.Rmd",
  output_dir = "/home/ocdm0351/DPhil/R_Data/htmls",
  quiet = FALSE
)'

echo "All jobs completed."

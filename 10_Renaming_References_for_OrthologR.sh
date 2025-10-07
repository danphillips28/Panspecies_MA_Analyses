#!/bin/bash

#SBATCH --job-name=10_Renaming_References_for_OrthologR
#SBATCH --array=1-7%2 # Should be 1 to the number of files (rows) in the summary file
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A_%a.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A_%a.err

# Paths
SOURCE_FILE="/home/ocdm0351/DPhil/liftingOver/liftOver_Source_File.txt"
SNPEFF_DATA_DIR="/home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-2/data"
LOG_DIR="/home/ocdm0351/DPhil/logs"

# Extract the line corresponding to this array task (skip header with +1)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SOURCE_FILE")

# Parse the row
IFS=$'\t' read -r file species start_build end_build chain_file_url chain_file_path end_build_fasta_url end_build_fasta_path end_build_gtf_url end_build_gtf_path end_build_cds_url end_build_cds_path end_build_aa_url end_build_aa_path build_check <<< "$LINE"

# Create species build directories
SPECIES_DIR="$SNPEFF_DATA_DIR/$species"
END_BUILD_DIR="$SPECIES_DIR/$end_build"

# Renaming cds
echo "Subsetting cds.fa to canonical isoforms..."
grep -f "$END_BUILD_DIR/cds.fa" file1.tabular > "$END_BUILD_DIR/cds_$end_build.fa"
#mv "$END_BUILD_DIR/cds.fa" "$END_BUILD_DIR/cds_$end_build.fa"

echo "Done with $end_build"

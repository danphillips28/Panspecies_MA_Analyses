#!/bin/bash
#SBATCH --job-name=4_building_snpEff_databases.sh
#SBATCH --output=/home/ocdm0351/DPhil/logs/4_snpeff_building_database_%A_%a.out
#SBATCH --error=/home/ocdm0351/DPhil/logs/4_snpeff_building_database_%A_%a.err
#SBATCH --array=1-8%1
#SBATCH --mem=200G

# Load snpEff environment
module load Anaconda3
source activate snpEff

SOURCE_FILE="/home/ocdm0351/DPhil/liftingOver/liftOver_Source_File.txt"
SNP_EFF_CONFIG="/home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/snpEff.config"
DATA_ROOT="/home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/data"

# Extract the line for this array task (skip header)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SOURCE_FILE)

# Parse the row using tab as delimiter
IFS=$'\t' read -r file species start_build end_build chain_file_url chain_file_path \
    end_build_fasta_url end_build_fasta_path end_build_gtf_url end_build_gtf_path \
    end_build_cds_url end_build_cds_path end_build_aa_url end_build_aa_path genome_codon_table mito_codon_table mito_chr build_check <<< "$LINE"

DATA_DIR="${DATA_ROOT}/${species}/${end_build}"
GENOME_NAME="${species}/${end_build}"

# Appending genome information to config file if not present
if ! grep -q "^${GENOME_NAME}\.genome" "$SNP_EFF_CONFIG"; then
    echo "Adding $GENOME_NAME to snpEff.config"
    {
        echo "${GENOME_NAME}.genome : ${species}"
        echo "${GENOME_NAME}.codonTable : ${genome_codon_table}"
        echo "${GENOME_NAME}.${mito_chr}.codonTable : ${mito_codon_table}"
    } >> "$SNP_EFF_CONFIG"
fi

# Build database
cd /home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/
snpEffDir="/home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1"
# Calling snpEff as below, using java, means we can use -Xmx to specify how much memory is used.
    java -Xmx160g -jar "$snpEffDir/snpEff.jar" build -gtf22 $build_check -v -c "$SNP_EFF_CONFIG" "$GENOME_NAME"

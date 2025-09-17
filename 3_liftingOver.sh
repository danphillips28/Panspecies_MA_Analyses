#!/bin/bash

#SBATCH --job-name=3_lifting_Over
#SBATCH --array=1-7%2 # Should be 1 to the number of files (rows) in the summary file
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A_%a.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A_%a.err

module load picard/3.0.0-Java-17
module load Anaconda3
source activate snpEff

# Paths
SOURCE_FILE="/home/ocdm0351/DPhil/liftingOver/liftOver_Source_File.txt"
SNPEFF_DATA_DIR="/home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/data"
CHAIN_DIR="/home/ocdm0351/DPhil/liftingOver/chainFiles"
OUTPUT_DIR="/home/ocdm0351/DPhil/LIFTED_VCF_DATA"
LOG_DIR="/home/ocdm0351/DPhil/logs"

# Make sure directories exist
mkdir -p "$SNPEFF_DATA_DIR" "$CHAIN_DIR" "$OUTPUT_DIR" "$LOG_DIR"

# Extract the line corresponding to this array task (skip header with +1)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SOURCE_FILE")

# Parse the row
IFS=$'\t' read -r file species start_build end_build chain_file_url chain_file_path end_build_fasta_url end_build_fasta_path end_build_gtf_url end_build_gtf_path end_build_cds_url end_build_cds_path end_build_aa_url end_build_aa_path build_check <<< "$LINE"

# Redirect stdout/stderr to log files named after input file
BASENAME="${file%.vcf}"
exec >"$LOG_DIR/${BASENAME}_lifting.log" 2>"$LOG_DIR/${BASENAME}_lifting.err"

echo "Processing: $file (species=$species, start=$start_build, end=$end_build)"

# Create species build directories
SPECIES_DIR="$SNPEFF_DATA_DIR/$species"
END_BUILD_DIR="$SPECIES_DIR/$end_build"
mkdir -p "$END_BUILD_DIR"

# Download references
echo "Downloading references..."
if [ -n "$end_build_fasta_url" ]; then wget -q -O "$END_BUILD_DIR/sequences.fa.gz" "$end_build_fasta_url" && gunzip -f "$END_BUILD_DIR/sequences.fa.gz"; fi
if [ -n "$end_build_gtf_url" ]; then wget -q -O "$END_BUILD_DIR/genes.gtf.gz" "$end_build_gtf_url" && gunzip -f "$END_BUILD_DIR/genes.gtf.gz"; fi
if [ -n "$end_build_cds_url" ]; then wget -q -O "$END_BUILD_DIR/cds.fa.gz" "$end_build_cds_url" && gunzip -f "$END_BUILD_DIR/cds.fa.gz"; fi
if [ -n "$end_build_aa_url" ]; then wget -q -O "$END_BUILD_DIR/protein.fa.gz" "$end_build_aa_url" && gunzip -f "$END_BUILD_DIR/protein.fa.gz"; fi

REFERENCE="$END_BUILD_DIR/sequences.fa"

# Make dictionary file if not exists
DICT_FILE="${REFERENCE%.fa}.dict"
if [ ! -f "$DICT_FILE" ]; then
    echo "Creating sequence dictionary..."
    java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R="$REFERENCE" O="$DICT_FILE"
fi

# Liftover or copy
INPUT_VCF="/home/ocdm0351/DPhil/Groomed_VCF_DATA/$file"
OUTPUT_VCF="$OUTPUT_DIR/${BASENAME}_${end_build}.vcf"

if [ "$start_build" != "$end_build" ]; then
    echo "Liftover required for $file ($start_build -> $end_build)"
    wget -q -P "$CHAIN_DIR/$species" "$chain_file_url"
    gunzip -f "$CHAIN_DIR/$species/$chain_file_path.gz"
    sed -i 's/chr//g' "$CHAIN_DIR/$species/$chain_file_path"
    CHAIN_FILE="$CHAIN_DIR/$species/$chain_file_path"
    REJECT_VCF="$OUTPUT_DIR/${BASENAME}_rejected.vcf"

    java -jar $EBROOTPICARD/picard.jar LiftoverVcf \
        I="$INPUT_VCF" \
        O="$OUTPUT_VCF" \
        CHAIN="$CHAIN_FILE" \
        REJECT="$REJECT_VCF" \
        R="$REFERENCE" \
        WARN_ON_MISSING_CONTIG=true
else
    echo "No liftover required, copying instead"
    cp "$INPUT_VCF" "$OUTPUT_VCF"
fi

echo "Done with $file"

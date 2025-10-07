#!/bin/bash

#SBATCH --job-name=9_Creating_Sister_Species_Table
#SBATCH --array=1-7 # Should be 1 to the number of files (rows) in the summary file
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A_%a.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A_%a.err

# First Narrow Down the liftOver_Source_File
#cat ../liftingOver/liftOver_Source_File.txt | cut -f 1-2,4,12 > ../R_Data/Sister_Species_Table.txt 

# Then Manually Add in the name, CDS URL, and CDS file for each sister species
# nano blah blah

module load Anaconda3
source activate snpEff

# Now for each row in the table, download and unzip the cds file
# Paths
SOURCE_FILE="/home/ocdm0351/DPhil/R_Data/Sister_Species_Table.txt"
SNPEFF_DATA_DIR="/home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-2/data"
LOG_DIR="/home/ocdm0351/DPhil/logs"

# Extract the line corresponding to this array task (skip header with +1)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$SOURCE_FILE")

# Parse the row
IFS=$'\t' read -r file species end_build end_build_cds_path sister_species sister_cds_url sister_cds_path <<< "$LINE"

# Redirect stdout/stderr to log files named after input file
BASENAME="$species_and_$sister_species"
exec >"$LOG_DIR/${BASENAME}_download.log" 2>"$LOG_DIR/${BASENAME}_download.err"

echo "Processing: (species=$species, sister_species=$sister_species)"

# Create species build directories
SISTER_SPECIES_DIR="$SNPEFF_DATA_DIR/$sister_species"
mkdir -p "$SISTER_SPECIES_DIR"

# Download references
echo "Downloading references..."
wget -q -O "$SISTER_SPECIES_DIR/cds.fa.gz" "$sister_cds_url" && gunzip -f "$SISTER_SPECIES_DIR/cds.fa"

echo "Done with sister species of $species: $sister_species"

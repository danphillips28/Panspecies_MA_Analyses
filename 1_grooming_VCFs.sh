#!/bin/bash
#SBATCH --job-name=GroomVCF   
#SBATCH --output=/home/ocdm0351/DPhil/logs/1_grooming_VCFs_%j.log  # where %j will be replaced by the job ID
#SBATCH --error=/home/ocdm0351/DPhil/logs/1_grooming_VCFs_%j.err   # optional: separate error file

# Input and output directories
INPUT_DIR="/home/ocdm0351/DPhil/VCF_DATA"
OUTPUT_DIR="/home/ocdm0351/DPhil/Groomed_VCF_DATA"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over all .vcf files in the input directory
for file in "$INPUT_DIR"/*.vcf; do
    # Extract the filename without path
    filename=$(basename "$file")
    # Define the output filename
    output_file="$OUTPUT_DIR/${filename%.vcf}_Groomed.vcf"

    # Check if the file already has a ##fileformat line
    if grep -q "^##fileformat" "$file"; then
        # If fileformat exists, just process the filters
	awk 'BEGIN{OFS="\t"} /^#/ {print; next} !($4==$5) && !(($4~/[ACGT]/) && ($5=="."))' "$file" > "$output_file"
    else
        # If fileformat is missing, add it at the top, then process the filters
        {
            echo "##fileformat=VCFv4.2"
            awk 'BEGIN{OFS="\t"} !($4==$5) && !(($4~/[ACGT]/) && ($5=="."))' "$file"
        } > "$output_file"
    fi

    echo "Processed $filename -> $(basename "$output_file")"
done

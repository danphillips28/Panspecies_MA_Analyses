#!/bin/bash
#SBATCH --job-name=2_GroomVCF
#SBATCH --output=/home/ocdm0351/DPhil/logs/2_grooming_VCFs_%j.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/2_grooming_VCFs_%j.err

INPUT_DIR="/home/ocdm0351/DPhil/VCF_DATA"
OUTPUT_DIR="/home/ocdm0351/DPhil/Groomed_VCF_DATA"
mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*.vcf; do
    filename=$(basename "$file")
    output_file="$OUTPUT_DIR/${filename%.vcf}_Groomed.vcf"

    # If the first line is not fileformat, add it to the top
    firstline=$(head -n 1 "$file")
    if [[ ! $firstline =~ ^##fileformat ]]; then
        echo "##fileformat=VCFv4.2" > "$output_file"
    else
        : > "$output_file"  # clear any previous contents
    fi

    # Then append the (cleaned) original content
    sed 's/[[:space:]]\+/\t/g' "$file" \
      | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {gsub("-", "", $4); gsub("-", "", $5); print}' \
      | awk 'BEGIN{OFS="\t"} !($4==$5)' \
      | awk 'BEGIN{OFS="\t"} !((length($4)==1) && ($4~/[ACGT]/) && ($5=="."))' \
      >> "$output_file"

    echo "Processed $filename -> $(basename "$output_file")"
done

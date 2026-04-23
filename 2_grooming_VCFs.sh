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

    # Ensure file starts with a fileformat header
    firstline=$(head -n 1 "$file")
    if [[ ! $firstline =~ ^##fileformat ]]; then
        echo "##fileformat=VCFv4.2" > "$output_file"
    else
        : > "$output_file"
    fi

    # Keep only SNPs:
    # - REF and ALT are one base long
    # - both are A/C/G/T
    # - ALT is not multiallelic
    # - REF != ALT
    awk 'BEGIN{OFS="\t"}
        /^#/ {print; next}
        {
            split($5, alts, ",")
            if (length($4)==1 && length($5)==1 &&
                $4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ &&
                $4 != $5) {
                print
            }
        }' "$file" >> "$output_file"

    echo "Processed $filename -> $(basename "$output_file")"
done

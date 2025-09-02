#!/bin/bash
#SBATCH --job-name=5_snpEffSift.sh
#SBATCH --output=/home/ocdm0351/DPhil/logs/5_snpEffSift_%A_%a.out
#SBATCH --error=/home/ocdm0351/DPhil/logs/5_snpEffSift_%A_%a.err
#SBATCH --array=1-7

module load Anaconda3
source activate snpEff

SOURCE_FILE="/home/ocdm0351/DPhil/liftingOver/liftOver_Source_File.txt"
SNP_EFF_CONFIG="/home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1/snpEff.config"
INPUT_DIREC="/home/ocdm0351/DPhil/LIFTED_VCF_DATA"
OUTPUT_DIREC="/home/ocdm0351/DPhil/snpEff/Annotated_VCFs"

# Extract the line for this array task (skip header)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SOURCE_FILE)

# Parse the row using tab as delimiter
IFS=$'\t' read -r file species start_build end_build chain_file_url chain_file_path \
    end_build_fasta_url end_build_fasta_path end_build_gtf_url end_build_gtf_path \
    end_build_cds_url end_build_cds_path end_build_aa_url end_build_aa_path \
    genome_codon_table mito_codon_table mito_chr build_check <<< "$LINE"


# Modifying the file name to represent the new post-liftover names
lifted_file="${file%.vcf}_$end_build.vcf"
GENOME_NAME="${species}/${end_build}"

# Run snpEff
snpEffDir="/home/ocdm0351/.conda/envs/snpEff/share/snpeff-5.2-1"
# Calling snpEff as below, using java, means we can use -Xmx to specify how much memory is used.
java -Xmx160g -jar $snpEffDir/snpEff.jar eff -noStats -v -d -c $SNP_EFF_CONFIG -canon \
		 $GENOME_NAME $INPUT_DIREC/$lifted_file  > "$OUTPUT_DIREC/Annotated_$lifted_file"

cat "$OUTPUT_DIREC/Annotated_$lifted_file" | /home/ocdm0351/DPhil/snpEff/vcfEffOnePerLine.pl | SnpSift extractFields - CHROM POS REF ALT FILTER "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].BIOTYPE" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" > "$OUTPUT_DIREC/OnePerLine_$lifted_file"

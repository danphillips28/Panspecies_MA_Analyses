#!/bin/bash

#SBATCH --job-name=17_Downloading_Dmelanogaster_Epigenetics_Covariates_Data
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A.err

# Load UCSC liftOver (conda environment)
module load Anaconda3
source activate UCSC_liftOver

input_file="/home/ocdm0351/DPhil/scripts/Dmelanogaster_Epigenetic_Data_Source_File.tsv"
species=$(head -1 "$input_file" | cut -f2 | sed 's/_URL//')

DATA_DIR="/home/ocdm0351/DPhil/R_Data"
CHAIN_DIR="/home/ocdm0351/DPhil/liftingOver/chain_files"
mkdir -p "$DATA_DIR" "$CHAIN_DIR/$species"

tail -n +2 "$input_file" | while IFS=$'\t' read -r modification url start_build end_build chain_file chain_file_path; do
    [[ -z "$url" ]] && continue

    outfile="$DATA_DIR/${species}_${modification}.bedgraph"
    lifted_file="${outfile%.bedgraph}_lifted.bedgraph"

    echo "Downloading and decompressing $outfile..."
    wget -q -O - "$url" | gunzip -c > "$outfile"

    # --- Simple cleanup: remove first line and replace spaces with tabs ---
    echo "Cleaning $outfile (remove header, convert spaces to tabs)..."
    sed '1d; s/ /\t/g' "$outfile" > "${outfile}.tmp" && mv "${outfile}.tmp" "$outfile"

    # --- Liftover step ---
    if [ "$start_build" != "$end_build" ]; then
        echo "Liftover required for $outfile ($start_build -> $end_build)"

        chain_name=$(basename "$chain_file_path")
        chain_path="$CHAIN_DIR/$species/$chain_name"

        # Download and decompress chain file if needed
        if [[ ! -f "${chain_path%.gz}" ]]; then
            echo "Downloading chain file..."
            wget -q -P "$CHAIN_DIR/$species" "$chain_file"
            if [[ "$chain_path" == *.gz ]]; then
                gunzip -f "$chain_path"
                chain_path="${chain_path%.gz}"
            fi
        fi

        unmapped_file="${outfile%.bedgraph}_unmapped.txt"

        echo "Running UCSC liftOver: liftOver "$outfile" "$chain_path" "$lifted_file" "$unmapped_file""
        liftOver "$outfile" "$chain_path" "$lifted_file" "$unmapped_file"

        if [[ -f "$unmapped_file" ]]; then
            mv -f "$unmapped_file" /home/ocdm0351/DPhil/logs/
        fi

        if [[ $? -eq 0 && -s "$lifted_file" ]]; then
            mv "$lifted_file" "$outfile"
            echo "Liftover completed successfully for $modification."
        else
            echo "Liftover failed or produced no lifted lines for $modification."
            rm -f "$lifted_file"
        fi
    else
        echo "No liftover required, keeping original file."
    fi

    echo "Done with $outfile"
    echo
done

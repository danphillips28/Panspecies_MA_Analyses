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

    tmp_download="$(mktemp)"
    tmp_unpacked="$(mktemp)"

    echo "Downloading $outfile..."
    if ! wget -q -O "$tmp_download" "$url"; then
        echo "Download failed for $modification"
        rm -f "$tmp_download" "$tmp_unpacked"
        continue
    fi

    echo "Inspecting compression for $outfile..."
    if file -b "$tmp_download" | grep -qi gzip; then
        gunzip -c "$tmp_download" > "$tmp_unpacked"
    else
        cp "$tmp_download" "$tmp_unpacked"
    fi

    echo "Cleaning and normalizing $outfile..."
    if [[ "$modification" == "replication_timing" ]]; then
        # Special case:
        # source file format is: chr, start, value
        # convert to: chr, start, end, value
        awk 'BEGIN { OFS="\t" }
            /^#/ { next }
            /^track/ { next }
            /^browser/ { next }
            $1 == "chr" && $2 == "start" { next }
            NF == 3 {
                print $1, $2, $2 + 1, $3
                next
            }
            NF >= 4 {
                print $1, $2, $3, $4
                next
            }' "$tmp_unpacked" > "$outfile"
    else
        # Standard files: keep the first four columns
        awk 'BEGIN { OFS="\t" }
            /^#/ { next }
            /^track/ { next }
            /^browser/ { next }
            $1 == "chr" && $2 == "start" { next }
            NF >= 4 {
                print $1, $2, $3, $4
                next
            }' "$tmp_unpacked" > "$outfile"
    fi

    rm -f "$tmp_download" "$tmp_unpacked"

    # --- Liftover step ---
    if [ "$start_build" != "$end_build" ]; then
        echo "Liftover required for $outfile ($start_build -> $end_build)"

        chain_name=$(basename "$chain_file_path")
        chain_path="$CHAIN_DIR/$species/$chain_name"

        # Download chain file if needed
        if [[ ! -f "$chain_path" && ! -f "${chain_path%.gz}" ]]; then
            echo "Downloading chain file..."
            if ! wget -q -O "$CHAIN_DIR/$species/$chain_name" "$chain_file"; then
                echo "Failed to download chain file for $modification"
                continue
            fi
        fi

        # Decompress chain file if needed
        if [[ -f "$chain_path" && "$chain_path" == *.gz ]]; then
            gunzip -f "$chain_path"
            chain_path="${chain_path%.gz}"
        elif [[ ! -f "$chain_path" && -f "${chain_path%.gz}" ]]; then
            chain_path="${chain_path%.gz}"
        fi

        unmapped_file="${outfile%.bedgraph}_unmapped.txt"

        echo "Running UCSC liftOver: liftOver \"$outfile\" \"$chain_path\" \"$lifted_file\" \"$unmapped_file\""
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

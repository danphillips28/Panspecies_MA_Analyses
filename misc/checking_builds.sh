#!/bin/bash

module load Anaconda3
source activate conda
source activate snpEff

# List your species here (exact names as in snpEff databases)
species_list=(
    Drosophila_melanogaster
    Schizosaccharomyces_pombe
    Neurospora_crassa
    Caenorhabditis_elegans
)

for sp in "${species_list[@]}"; do
    echo "=== $sp ==="
    java -jar snpEff.jar download -v "$sp" 2>&1 \
        | grep -E "ftp|http" \
        | grep -i -E "ensembl|ncbi|ucsc" \
        | head -n 1
    echo
done

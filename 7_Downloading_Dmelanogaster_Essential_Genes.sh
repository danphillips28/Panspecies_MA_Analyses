#!/bin/bash
set -euo pipefail

PHENO_URL="https://s3ftp.flybase.org/releases/current/precomputed_files/alleles/genotype_phenotype_data_fb_2026_01.tsv.gz"
MAP_URL="https://s3ftp.flybase.org/releases/current/precomputed_files/alleles/fbal_to_fbgn_fb_2026_01.tsv.gz"

OUT="/home/ocdm0351/DPhil/R_Data/Dmelanogaster_Essential_Genes.tsv"

wget -q -O phen.tsv.gz "$PHENO_URL"
wget -q -O map.tsv.gz  "$MAP_URL"

# Step 1: get alleles with lethal/inviable/die
zcat phen.tsv.gz | grep -v '^#' | awk -F'\t' '
  tolower($3) ~ /(lethal|inviable|die)/ && tolower($3) !~ /\bviable\b/ {print $2}
' | sort -u > alleles.txt

# Step 2: map alleles -> FBgn (robustly grab FBgn column)
zcat map.tsv.gz | awk -F'\t' '
  NR==FNR {a[$1]; next}
  {
    for(i=1;i<=NF;i++){
      if($i ~ /^FBgn/ && $1 in a){
        print $i
      }
    }
  }
' alleles.txt - | sort -u > "$OUT"

echo "Done: $OUT"

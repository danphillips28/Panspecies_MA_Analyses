#!/bin/bash

#SBATCH --job-name=15_Downloading_Ppacificus_Epigenetics_Covariates_Data
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A.err

module load all/deepTools/3.5.2-foss-2022a
module load all/BEDTools/2.31.0-GCC-12.3.0

set -euo pipefail
IFS=$'\n\t'   # keep newline+tab for most parsing

TSV="Ppacificus_Epigenetic_Data_Source_File.tsv"
OUTDIR="/home/ocdm0351/DPhil/R_Data"
BIN_SIZE=25
THREADS=4

mkdir -p "$OUTDIR"

# group URLs by mark into MARK<TAB>url1|url2|...
awk -F'\t' 'NR>1 {g=$1; url=$2; gsub(/ /,"_",g); a[g]=a[g]?a[g]"|"url:url} END{for(x in a) print x "\t" a[x]}' "$TSV" > "$OUTDIR/_marks.tsv"

while IFS=$'\t' read -r MARK URLS; do
  echo "Processing $MARK"

  # MINIMAL CHANGE: skip if final bedgraph already exists
  final="${OUTDIR}/P_pacificus_${MARK}.bedgraph"
  if [ -f "$final" ]; then
    echo "Skipping $MARK (already processed): $final"
    continue
  fi
  # End minimal change

  MDIR="$OUTDIR/$MARK"; mkdir -p "$MDIR"

  # split the pipe-separated URL list into an array
  IFS='|' read -r -a url_array <<< "$URLS"

  bgs=(); i=0
  for url in "${url_array[@]}"; do
    i=$((i+1))
    bamfn="${MDIR}/$(basename "$url")"
    [[ "$bamfn" != *.bam ]] && bamfn="${MDIR}/${MARK}_rep${i}.bam"

    echo "  downloading $url -> $bamfn"
    wget -c -O "$bamfn" "$url"
    # assume .bai available at URL.bai:
    wget -c -O "${bamfn}.bai" "${url}.bai" || true

    bg="${bamfn%.bam}.bedgraph"
    bamCoverage -b "$bamfn" -o "$bg" --outFileFormat bedgraph --normalizeUsing CPM --binSize "$BIN_SIZE" --numberOfProcessors "$THREADS"

    bgs+=("$bg")
  done

  if [ "${#bgs[@]}" -eq 1 ]; then
    cp -f "${bgs[0]}" "$final"
  else
    bedtools unionbedg -i "${bgs[@]}" > "${MDIR}/union.tab"
    awk -v OFS="\t" '{ s=0;n=0; for(i=4;i<=NF;i++){ if($i!="" && $i!="nan"){ s+=($i); n++ } } mean=(n? s/n:0); print $1,$2,$3,mean }' "${MDIR}/union.tab" > "$final"
  fi
  echo "Wrote $final"

	rm -rf "$MDIR"

done < "$OUTDIR/_marks.tsv"

echo "All done. Bedgraphs in: $OUTDIR"

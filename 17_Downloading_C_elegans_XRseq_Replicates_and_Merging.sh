#!/usr/bin/env bash

#SBATCH --job-name=17_Downloading_C_elegans_XRseq_Replicates_and_Merging.sh
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A.err
#SBATCH --time=02:00:00
#SBATCH --partition=himem-gen
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

module load BEDTools

set -euo pipefail

OUTDIR="/home/ocdm0351/DPhil/R_Data/"
MANIFEST="/home/ocdm0351/DPhil/scripts/Cel_XRseq_Manifest.tsv"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

if ! command -v bedtools >/dev/null 2>&1; then
  echo "Error: bedtools not found in PATH" >&2
  exit 1
fi

declare -A FILELIST
declare -A COUNT
declare -A SEEN
SAMPLE_ORDER=()

while IFS=$'\t' read -r SAMPLE STRAND URL; do
  [[ -z "${SAMPLE:-}" ]] && continue
  [[ "${SAMPLE}" =~ ^# ]] && continue
  [[ "${SAMPLE}" == "sample_name" ]] && continue

  if [[ "${STRAND}" != "minus" && "${STRAND}" != "plus" ]]; then
    echo "Error: strand must be minus or plus for sample ${SAMPLE}" >&2
    exit 1
  fi

  key="${SAMPLE}|${STRAND}"

  if [[ -z "${SEEN["$SAMPLE"]+x}" ]]; then
    SAMPLE_ORDER+=("${SAMPLE}")
    SEEN["$SAMPLE"]=1
  fi

  COUNT["$key"]=$(( ${COUNT["$key"]:-0} + 1 ))
  rep="${COUNT["$key"]}"

  SAFE_SAMPLE=$(printf '%s' "${SAMPLE}" | tr -cs 'A-Za-z0-9._-' '_')
  gzfile="${SAFE_SAMPLE}_${STRAND}_rep${rep}.bedgraph.gz"
  bedfile="${SAFE_SAMPLE}_${STRAND}_rep${rep}.bedgraph"

  echo "Downloading ${SAMPLE} ${STRAND} rep${rep}"
  wget -O "${gzfile}" "${URL}"

  echo "Unzipping ${gzfile}"
  gunzip -c "${gzfile}" \
    | awk 'BEGIN{OFS="\t"} !/^track/ && !/^browser/ {print}' \
    | sort -k1,1 -k2,2n \
    > "${bedfile}"

  if [[ -z "${FILELIST["$key"]+x}" ]]; then
    FILELIST["$key"]="${bedfile}"
  else
    FILELIST["$key"]+=$'\n'"${bedfile}"
  fi
done < "${MANIFEST}"

for SAMPLE in "${SAMPLE_ORDER[@]}"; do
  SAFE_SAMPLE=$(printf '%s' "${SAMPLE}" | tr -cs 'A-Za-z0-9._-' '_')

  for STRAND in minus plus; do
    key="${SAMPLE}|${STRAND}"
    [[ -z "${FILELIST["$key"]+x}" ]] && continue

    mapfile -t INPUTS <<< "${FILELIST["$key"]}"
    # Keep the underscore before plus/minus so downstream scripts can match *_plus$ and *_minus$
    out_bed="${SAFE_SAMPLE}_${STRAND}.bedgraph"

    echo "Merging ${SAMPLE} ${STRAND} -> ${out_bed}"

    if [[ "${#INPUTS[@]}" -eq 1 ]]; then
      cp "${INPUTS[0]}" "${out_bed}"
    else
      bedtools unionbedg -i "${INPUTS[@]}" -filler 0 \
        | awk 'BEGIN{OFS="\t"} {sum=0; n=NF-3; for (i=4; i<=NF; i++) sum += $i; print $1, $2, $3, sum/n}' \
        > "${out_bed}"
    fi
  done
done

echo "Done!"

rm /home/ocdm0351/DPhil/R_Data/*C_elegans*rep*

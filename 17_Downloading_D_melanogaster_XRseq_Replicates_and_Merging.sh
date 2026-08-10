#!/usr/bin/env bash
#SBATCH --job-name=17_Downloading_D_melanogaster_XRseq_Replicates_and_Merging.sh
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A.err
#SBATCH --time=02:00:00
#SBATCH --partition=himem-gen
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail

module load deepTools/3.5.2-foss-2022a
module load Anaconda3
source activate UCSC_liftOver

OUTDIR="/home/ocdm0351/DPhil/R_Data/"
cd "${OUTDIR}"

MANIFEST="/home/ocdm0351/DPhil/scripts/Dmel_XRseq_Manifest.tsv"

declare -A MINUS_FILES
declare -A PLUS_FILES
declare -A SEEN_SAMPLES

while IFS=$'\t' read -r SAMPLE STRAND URL; do
  [[ -z "${SAMPLE}" ]] && continue
  [[ "${SAMPLE}" =~ ^# ]] && continue

  SEEN_SAMPLES["${SAMPLE}"]=1

  FILENAME="${SAMPLE}_${STRAND}_$(basename "${URL}")"
  wget -O "${FILENAME}" "${URL}"

  if [[ "${STRAND}" == "minus" ]]; then
    if [[ -z "${MINUS_FILES[${SAMPLE}]+x}" ]]; then
      MINUS_FILES["${SAMPLE}"]="${FILENAME}"
    else
      MINUS_FILES["${SAMPLE}"]+=" ${FILENAME}"
    fi
  elif [[ "${STRAND}" == "plus" ]]; then
    if [[ -z "${PLUS_FILES[${SAMPLE}]+x}" ]]; then
      PLUS_FILES["${SAMPLE}"]="${FILENAME}"
    else
      PLUS_FILES["${SAMPLE}"]+=" ${FILENAME}"
    fi
  else
    echo "Unknown strand '${STRAND}' for sample '${SAMPLE}'" >&2
    exit 1
  fi
done < "${MANIFEST}"

for SAMPLE in "${!SEEN_SAMPLES[@]}"; do
  echo "Processing ${SAMPLE}..."

  MINUS_LIST=${MINUS_FILES["${SAMPLE}"]:-}
  PLUS_LIST=${PLUS_FILES["${SAMPLE}"]:-}

  if [[ -n "${MINUS_LIST}" ]]; then
    set -- ${MINUS_LIST}
    if [[ $# -eq 1 ]]; then
      cp "$1" "${SAMPLE}_XRminus.bw"
    else
      bigwigCompare -b1 "$1" -b2 "$2" --operation mean -o "${SAMPLE}_XRminus.bw" --binSize 25
      shift 2
      while [[ $# -gt 0 ]]; do
        TMP="${SAMPLE}_XRminus_tmp.bw"
        bigwigCompare -b1 "${SAMPLE}_XRminus.bw" -b2 "$1" --operation mean -o "${TMP}" --binSize 25
        mv "${TMP}" "${SAMPLE}_XRminus.bw"
        shift
      done
    fi
    bigWigToBedGraph "${SAMPLE}_XRminus.bw" "${SAMPLE}_XR_minus.bedgraph"
  fi

  if [[ -n "${PLUS_LIST}" ]]; then
    set -- ${PLUS_LIST}
    if [[ $# -eq 1 ]]; then
      cp "$1" "${SAMPLE}_XRplus.bw"
    else
      bigwigCompare -b1 "$1" -b2 "$2" --operation mean -o "${SAMPLE}_XRplus.bw" --binSize 25
      shift 2
      while [[ $# -gt 0 ]]; do
        TMP="${SAMPLE}_XRplus_tmp.bw"
        bigwigCompare -b1 "${SAMPLE}_XRplus.bw" -b2 "$1" --operation mean -o "${TMP}" --binSize 25
        mv "${TMP}" "${SAMPLE}_XRplus.bw"
        shift
      done
    fi
    bigWigToBedGraph "${SAMPLE}_XRplus.bw" "${SAMPLE}_XR_plus.bedgraph"
  fi
done

echo "Done!"

rm /home/ocdm0351/DPhil/R_Data/*D_melanogaster*rep*
rm /home/ocdm0351/DPhil/R_Data/*.bw*

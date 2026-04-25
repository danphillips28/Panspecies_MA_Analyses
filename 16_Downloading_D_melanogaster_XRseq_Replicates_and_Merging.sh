#!/usr/bin/env bash
set -euo pipefail

module load deepTools/3.5.2-foss-2022a
module load Anaconda3
source activate UCSC_liftOver

# Config
OUTDIR="/home/ocdm0351/DPhil/R_Data/"
REP1minus_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5600nnn/GSM5600478/suppl/GSM5600478%5FS2CPD1h1%5FACAGTG%5FS9%5FL001%5FR1%5F001%5Fminus.bw"
REP2minus_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5661nnn/GSM5661670/suppl/GSM5661670%5FS2CPD1h2%5FGCCAAT%5FS10%5FL001%5FR1%5F001%5Fminus.bw"
REP1plus_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5600nnn/GSM5600478/suppl/GSM5600478%5FS2CPD1h1%5FACAGTG%5FS9%5FL001%5FR1%5F001%5Fplus.bw"
REP2plus_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5661nnn/GSM5661670/suppl/GSM5661670%5FS2CPD1h2%5FGCCAAT%5FS10%5FL001%5FR1%5F001%5Fplus.bw"

cd /home/ocdm0351/DPhil/R_Data/

# Download
echo "Downloading replicates..."
wget -O Dmelanogaster_XR1minus.bw "${REP1minus_URL}"
wget -O Dmelanogaster_XR2minus.bw "${REP2minus_URL}"
wget -O Dmelanogaster_XR1plus.bw "${REP1plus_URL}"
wget -O Dmelanogaster_XR2plus.bw "${REP2plus_URL}"

# Merge (average)
echo "Averaging Minus Strand replicates into one consensus bigWig..."
bigwigCompare \
  -b1 /home/ocdm0351/DPhil/R_Data/Dmelanogaster_XR1minus.bw \
  -b2 /home/ocdm0351/DPhil/R_Data/Dmelanogaster_XR2minus.bw \
  --operation mean \
  -o Dmelanogaster_XRminus.bw \
  --binSize 25 

echo "Averaging Plus Strand replicates into one consensus bigWig..."
bigwigCompare \
  -b1 /home/ocdm0351/DPhil/R_Data/Dmelanogaster_XR1plus.bw \
  -b2 /home/ocdm0351/DPhil/R_Data/Dmelanogaster_XR2plus.bw \
  --operation mean \
  -o Dmelanogaster_XRplus.bw \
  --binSize 25 
echo "Done!"

# 1) bigWig -> bedGraph  (no chrom.sizes needed)
echo "Converting bigWig to bedGraph..."
bigWigToBedGraph Dmelanogaster_XRminus.bw D_melanogaster_XRminus.bedgraph
bigWigToBedGraph Dmelanogaster_XRplus.bw D_melanogaster_XRplus.bedgraph

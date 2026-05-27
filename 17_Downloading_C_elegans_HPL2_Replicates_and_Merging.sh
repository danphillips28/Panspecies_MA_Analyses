#!/usr/bin/env bash
set -euo pipefail

#SBATCH --job-name=17_Downloading_C_elegans_HPL2_Replicates_and_Merging
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A.err

module load deepTools/3.5.2-foss-2022a
module load Anaconda3
source activate UCSC_liftOver

# Config
OUTDIR="/home/ocdm0351/DPhil/R_Data/"
REP1_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2333nnn/GSM2333101/suppl/GSM2333101%5FHPL2%5FYA%5FChIPseq%5FBEADSmapq0%5Frep1.bw"
REP2_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2333nnn/GSM2333102/suppl/GSM2333102%5FHPL2%5FYA%5FChIPseq%5FBEADSmapq0%5Frep2.bw"

cd /home/ocdm0351/DPhil/R_Data/

# Download
echo "Downloading replicates..."
wget -O Celegans_HPL2_rep1.bw "${REP1_URL}"
wget -O Celegans_HPL2_rep2.bw "${REP2_URL}"

# QC: correlation
echo "Running replicate QC..."
multiBigwigSummary bins -b Celegans_HPL2_rep1.bw Celegans_HPL2_rep1.bw \
  --binSize 10000 \
  -o Celegans_HPL2_summary.npz \
  --outRawCounts Celegans_HPL2_rawCounts.tab

plotCorrelation -in Celegans_HPL2_summary.npz \
  --corMethod pearson \
  --whatToPlot heatmap \
  --plotNumbers \
  --plotTitle "Replicate correlation – HPL2_YA" \
  -o Celegans_HPL2_replicate_correlation.png \
  --labels rep1 rep2

#echo "Check replicate_correlation.png and rawCounts.tab to confirm high correlation."

# Merge (average)
echo "Averaging replicates into one consensus bigWig..."
bigwigCompare \
  -b1 /home/ocdm0351/DPhil/R_Data/Celegans_HPL2_rep1.bw \
  -b2 /home/ocdm0351/DPhil/R_Data/Celegans_HPL2_rep2.bw \
  --operation mean \
  -o Celegans_HPL2.bigwig \
  --binSize 25 

echo "Done!"

# 1) bigWig -> bedGraph  (no chrom.sizes needed)
echo "Converting bigWig to bedGraph..."
bigWigToBedGraph Celegans_HPL2.bigwig Celegans_HPL2_prelift.bedgraph

# 2) LiftOver  (no chrom.sizes needed)
echo "Downloading liftOver chain file"

echo "Running UCSC liftOver: liftOver on Celegans_HPL2.bigwig"
liftOver Celegans_HPL2_prelift.bedgraph \
        /home/ocdm0351/DPhil/liftingOver/chain_files/C_elegans/ce10ToCe11.over.chain \
        C_elegans_HPL2.bedgraph \
        /home/ocdm0351/DPhil/logs/unmapped_Celegans_HPL2.bedgraph

rm /home/ocdm0351/DPhil/R_Data/*prelift*

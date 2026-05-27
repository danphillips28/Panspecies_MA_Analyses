#!/usr/bin/env bash
set -euo pipefail

#SBATCH --job-name=17_Downloading_C_elegans_H3K9me2_Replicates_and_Merging
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A.err

module load deepTools/3.5.2-foss-2022a
module load Anaconda3
source activate UCSC_liftOver

# Config
OUTDIR="/home/ocdm0351/DPhil/R_Data/"
REP1_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2333nnn/GSM2333097/suppl/GSM2333097_H3K9me2_YA_ChIPseq_BEADSmapq0_rep1.bw"
REP2_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2333nnn/GSM2333098/suppl/GSM2333098_H3K9me2_YA_ChIPseq_BEADSmapq0_rep2.bw"

cd /home/ocdm0351/DPhil/R_Data/

# Download
echo "Downloading replicates..."
wget -O Celegans_H3K9me2_rep1.bw "${REP1_URL}"
wget -O Celegans_H3K9me2_rep2.bw "${REP2_URL}"

# QC: correlation
echo "Running replicate QC..."
multiBigwigSummary bins -b Celegans_H3K9me2_rep1.bw Celegans_H3K9me2_rep1.bw \
  --binSize 10000 \
  -o Celegans_H3K9me2_summary.npz \
  --outRawCounts Celegans_H3K9me2_rawCounts.tab

plotCorrelation -in Celegans_H3K9me2_summary.npz \
  --corMethod pearson \
  --whatToPlot heatmap \
  --plotNumbers \
  --plotTitle "Replicate correlation – H3K9me2_YA" \
  -o Celegans_H3K9me2_replicate_correlation.png \
  --labels rep1 rep2

#echo "Check replicate_correlation.png and rawCounts.tab to confirm high correlation."

# Merge (average)
echo "Averaging replicates into one consensus bigWig..."
bigwigCompare \
  -b1 /home/ocdm0351/DPhil/R_Data/Celegans_H3K9me2_rep1.bw \
  -b2 /home/ocdm0351/DPhil/R_Data/Celegans_H3K9me2_rep2.bw \
  --operation mean \
  -o Celegans_H3K9me2.bigwig \
  --binSize 25 

echo "Done!"

# 1) bigWig -> bedGraph  (no chrom.sizes needed)
echo "Converting bigWig to bedGraph..."
bigWigToBedGraph Celegans_H3K9me2.bigwig Celegans_H3K9me2_prelift.bedgraph

# 2) LiftOver  (no chrom.sizes needed)
echo "Downloading liftOver chain file"

wget -q -O - "https://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz" \
  | gunzip -c \
  > /home/ocdm0351/DPhil/liftingOver/chain_files/C_elegans/ce10ToCe11.over.chain

echo "Running UCSC liftOver: liftOver on Celegans_H3K9me2.bigwig"
liftOver Celegans_H3K9me2_prelift.bedgraph \
	/home/ocdm0351/DPhil/liftingOver/chain_files/C_elegans/ce10ToCe11.over.chain \
	C_elegans_H3K9me2.bedgraph \
	/home/ocdm0351/DPhil/logs/unmapped_Celegans_H3K9me2.bedgraph

rm /home/ocdm0351/DPhil/R_Data/*prelift*

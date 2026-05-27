#!/usr/bin/env bash


#SBATCH --job-name=17_Downloading_C_elegans_TimeSeries_Expression2
#SBATCH --output=/home/ocdm0351/DPhil/logs/%x_%A.log
#SBATCH --error=/home/ocdm0351/DPhil/logs/%x_%A.err


set -euo pipefail
trap 'echo "❌ Error on line $LINENO"' ERR

#wget -q -O - \
#https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50548/suppl/GSE50548%5FWhole%5Fembryo%5Finterval%5Ftimecourse.tab.gz \
#| gunzip -c > "/home/ocdm0351/DPhil/scripts/misc/GSE50548_Whole_embryo_interval_timecourse.tab"

# Cut just gene names to upload to wormbase gene sanitizer
#cat /home/ocdm0351/DPhil/scripts/misc/GSE50548_Whole_embryo_interval_timecourse.tab | cut -f 1 | sort | uniq \
#			 > /home/ocdm0351/DPhil/R_Data/GSE50548_Whole_embryo_interval_timecourse_gene_names.txt

# Used WormBase Gene Sanitizer to Create Mapping File and Upload to Genoa
grep -o '<tr>.*</tr>' /home/ocdm0351/DPhil/scripts/misc/gene_sanitizer.cgi.html \
	| sed -e 's/<\/t[dh]>/\t/g' -e 's/<[^>]*>//g' \
	| awk 'NF' > /home/ocdm0351/DPhil/scripts/misc/gene_sanitizer.cgi.tsv

# ==========================================================
# Replace gene names (first column) in GEO time series file
# with Ensembl/WBGene IDs from WormBase Gene Sanitizer output
# ==========================================================

# Input files
MAPPING="/home/ocdm0351/DPhil/scripts/misc/gene_sanitizer.cgi.tsv"
TIMESERIES="/home/ocdm0351/DPhil/scripts/misc/GSE50548_Whole_embryo_interval_timecourse.tab"
OUTFILE="/home/ocdm0351/DPhil/scripts/misc/GSE50548_WS235_converted.tab"
OUTDIR="/home/ocdm0351/DPhil/R_Data"

# Replace first column (gene name) with Ensembl/WBGene ID from mapping file
# Assumes: mapping file column 1 = old name, column 3 = new Ensembl/WBGene ID
awk -F'\t' -v OFS='\t' '
  NR==FNR { if (NF>=3) map[$1]=$3; next }
  NR==1 { print; next }  # keep header unchanged
  {
    $1 = (map[$1] != "") ? map[$1] : $1
    print
  }' "$MAPPING" "$TIMESERIES" > "$OUTFILE"

#echo "✅ Done. Output written to:"
#echo "   $OUTFILE"


echo "=== Splitting columns from $OUTFILE into $OUTDIR ==="

# Create files: column index -> desired suffix
# Column 1 = gene id, data columns start at 2
cut -f 1,2  "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_1C_Expression.tsv"
cut -f 1,3  "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_2C_Expression.tsv"
cut -f 1,4  "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_010_Expression.tsv"
cut -f 1,5  "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_020_Expression.tsv"
cut -f 1,6  "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_040_Expression.tsv"
cut -f 1,7  "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_050_Expression.tsv"
cut -f 1,8  "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_060_Expression.tsv"
cut -f 1,9  "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_070_Expression.tsv"
cut -f 1,10 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_100_Expression.tsv"
cut -f 1,11 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_120_Expression.tsv"
cut -f 1,12 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_130_Expression.tsv"
cut -f 1,13 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_150_Expression.tsv"
cut -f 1,14 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_160_Expression.tsv"
cut -f 1,15 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_170_Expression.tsv"
cut -f 1,16 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_180_Expression.tsv"
cut -f 1,17 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_190_Expression.tsv"
cut -f 1,18 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_200_Expression.tsv"
cut -f 1,19 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_210_Expression.tsv"
cut -f 1,20 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_220_Expression.tsv"
cut -f 1,21 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_230_Expression.tsv"
cut -f 1,22 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_260_Expression.tsv"
cut -f 1,23 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_270_Expression.tsv"
cut -f 1,24 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_280_Expression.tsv"
cut -f 1,25 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_290_Expression.tsv"
cut -f 1,26 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_300_Expression.tsv"
cut -f 1,27 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_320_Expression.tsv"
cut -f 1,28 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_330_Expression.tsv"
cut -f 1,29 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_340_Expression.tsv"
cut -f 1,30 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_350_Expression.tsv"
cut -f 1,31 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_360_Expression.tsv"
cut -f 1,32 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_370_Expression.tsv"
cut -f 1,33 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_380_Expression.tsv"
cut -f 1,34 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_390_Expression.tsv"
cut -f 1,35 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_400_Expression.tsv"
cut -f 1,36 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_410_Expression.tsv"
cut -f 1,37 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_420_Expression.tsv"
cut -f 1,38 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_430_Expression.tsv"
cut -f 1,38 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_Expression.tsv"
cut -f 1,39 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_450_Expression.tsv"
cut -f 1,40 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_490_Expression.tsv"
cut -f 1,41 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_510_Expression.tsv"
cut -f 1,42 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_550_Expression.tsv"
cut -f 1,43 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_620_Expression.tsv"
cut -f 1,44 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_630_Expression.tsv"
cut -f 1,45 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_640_Expression.tsv"
cut -f 1,46 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_700_Expression.tsv"
cut -f 1,47 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_750_Expression.tsv"
cut -f 1,48 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_760_Expression.tsv"
cut -f 1,49 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_800_Expression.tsv"
cut -f 1,50 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_810_Expression.tsv"
cut -f 1,51 "$OUTFILE" > "$OUTDIR/C_elegans_Embryonic_TimeSeries_830_Expression.tsv"

echo "All files written to $OUTDIR"

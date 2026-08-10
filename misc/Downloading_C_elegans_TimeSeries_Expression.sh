#!/bin/bash
set -euo pipefail

MISC="/home/ocdm0351/DPhil/scripts/misc"
OUT="${MISC}/Celegans_TimeSeries_Data_From_Boeck_etal.tsv"
MAP="${MISC}/c_elegans_clone_name_mapping.tsv"

echo "🧩 Downloading Ensembl Metazoa GTF (contains clone names)..."
wget -q -O - \
https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-59/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.59.gtf.gz \
| gunzip -c > "${MISC}/Caenorhabditis_elegans.WBcel235.59.gtf"

echo "✅ Building gene-name → WBGene mapping..."
awk '
  $3=="gene" {
    id=""; name=""; seq="";
    if(match($0,/gene_id "([^"]+)"/,a)) id=a[1];
    if(match($0,/gene_name "([^"]+)"/,b)) name=b[1];
    if(match($0,/sequence_name "([^"]+)"/,c)) seq=c[1];
    if(id && name) print name "\t" id;
    if(id && seq)  print seq  "\t" id;
  }
' "${MISC}/Caenorhabditis_elegans.WBcel235.59.gtf" | sort -u > "$MAP"

echo "✅ Mapping file created: $(wc -l < "$MAP") entries"

echo "📥 Downloading and remapping Boeck et al. expression dataset..."

wget -q -O - \
"https://genome.cshlp.org/content/suppl/2016/09/20/gr.202663.115.DC1/Supplemental_Table_S2.gz" \
| gunzip -c \
| awk 'BEGIN{FS="[ \t]+";OFS="\t"}{$1=$1;gsub(/[#"]/,"",$1);print}' \
| awk -v mapfile="$MAP" '
  BEGIN{
    FS=OFS="\t";
    while((getline < mapfile)>0){ map[$1]=$2 }
  }
  FNR==1 {print; next}
  { if($1 in map) $1=map[$1]; print }
' > "$OUT"

echo "🎉 Done!"
echo "Output file with Ensembl/WBGene IDs: $OUT"

# Extract the First Embryonic Timepoint
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,6 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_25_Expression.tsv"

# And so on..
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,8 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_55_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,10 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_85_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,12 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_235_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,14 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_265_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,16 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_325_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,18 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_355_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,20 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_385_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,22 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_415_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,24 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_445_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,26 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_475_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,28 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_505_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,30 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_535_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,32 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_565_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,34 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_595_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,36 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_625_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,38 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_655_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,40 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_685_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,42 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_715_Expression.tsv"
cat /home/ocdm0351/DPhil/scripts/misc/Celegans_TimeSeries_Data_From_Boeck_etal.tsv | cut -f 1,44 > "/home/ocdm0351/DPhil/R_Data/C_elegans_Embryonic_TimeSeries_745_Expression.tsv"

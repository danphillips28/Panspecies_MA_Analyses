#!/bin/bash

# Download the entire dataset
wget -q -O - "https://s3ftp.flybase.org/releases/current/precomputed_files/genes/gene_rpkm_report_fb_2025_04.tsv.gz" | gunzip -c > "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv"

# Extract the First Embryonic Timepoint
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em0-2hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_0_Expression.tsv"

# And so on..
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em2-4hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_2_Expression.tsv"
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em4-6hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_4_Expression.tsv"
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em6-8hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_6_Expression.tsv"
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em8-10hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_8_Expression.tsv"
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em10-12hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_10_Expression.tsv"
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em12-14hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_12_Expression.tsv"
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em14-16hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_14_Expression.tsv"
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em14-16hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_Expression.tsv"
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em16-18hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_16_Expression.tsv"
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em18-20hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_18_Expression.tsv"
awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em22-24hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" | cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_TimeSeries_22_Expression.tsv"

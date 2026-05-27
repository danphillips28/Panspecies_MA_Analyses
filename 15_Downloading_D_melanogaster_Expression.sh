#!/bin/bash 

wget -q -O - "https://s3ftp.flybase.org/releases/current/precomputed_files/genes/gene_rpkm_report_fb_2025_04.tsv.gz" \
		| gunzip -c > "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" 

awk -F'\t' 'NR==6 || $7 == "mE_mRNA_em14-16hr"' "/home/ocdm0351/DPhil/scripts/misc/Drosophila_melanogaster_modENCODE_RPKM_REPORT.tsv" \
		| cut -f2,8 > "/home/ocdm0351/DPhil/R_Data/D_melanogaster_Embryonic_Expression.tsv"

#!/bin/bash
wget -q -O - "https://www.encodeproject.org/files/ENCFF736CCO/@@download/ENCFF736CCO.tsv" \
  | cut -f1,7 \
  > /home/ocdm0351/DPhil/R_Data/H_sapiens_Embryonic_Expression.tsv

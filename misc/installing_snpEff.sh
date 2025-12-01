#!/bin/bash

# Load Anaconda3 module
module load Anaconda3
conda create -n snpEff
source activate snpEff
conda install mamba -c conda-forge
mamba install bioconda::snpeff bioconda::snpsift bioconda::seqtk	


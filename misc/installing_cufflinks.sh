#!/bin/bash

# Load Anaconda3 module
module load Anaconda3
conda create -n cufflinks
source activate cufflinks
conda install mamba -c conda-forge
mamba install bioconda/label/cf201901::cufflinks

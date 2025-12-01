#!/bin/bash

module load Anaconda3
conda create -n UCSC_liftOver
source activate UCSC_liftOver
conda install bioconda::ucsc-liftover

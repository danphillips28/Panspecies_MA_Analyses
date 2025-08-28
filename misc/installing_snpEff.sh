#!/bin/bash

# Load Anaconda3 module
module load Anaconda3
conda create -n snpEff python=3.11
conda activate snpEff
conda install mamba -c conda-forge

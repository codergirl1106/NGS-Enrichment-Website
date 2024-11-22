#!/bin/bash

# Update package list and install R and dependencies
sudo apt-get update && sudo apt-get install -y \
    r-base \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    build-essential

# Install BiocManager in R
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"

# Install required R packages
Rscript -e "BiocManager::install(c('Biostrings', 'Rfastp'))"

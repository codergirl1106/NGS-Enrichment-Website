#!/bin/bash

# Install BiocManager if not already installed
R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"

# Use BiocManager to install the necessary Bioconductor packages
R -e "BiocManager::install(c('Biostrings', 'Rfastp'))"

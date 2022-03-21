#!/bin/bash

# Create conda environment
mkdir envs
env=envs/default
conda create --prefix $env -y

# Activate environment and install packages
source activate $env

conda install -c conda-forge r-base=4.1.2 -y
conda install -c conda-forge r-seurat=4.1.0 -y
conda install -c bioconda bioconductor-rhdf5=2.38.0 -y

# Rscript -e 'install.packages("package", repos = "https://cloud.r-project.org/")'

# Create hidden file to indicate active environment
echo $env > .env

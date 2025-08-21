# syntax=docker/dockerfile:1

# FROM rocker/r-ver:4.4.1 

FROM satijalab/seurat:latest

# libglpk40 is needed for Seurat
RUN apt-get update && apt-get install -y libglpk40 && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -q -e "install.packages('Seurat')"
RUN R -q -e "install.packages(c('dplyr', 'ggplot2'))"
RUN R -q -e "install.packages('cowplot')"


# Get the preprocessing script and make it executable
COPY scripts/preprocess_seurat.R scripts/preprocess_seurat_utils.R /usr/local/bin/
RUN chmod +x /usr/local/bin/preprocess_seurat.R

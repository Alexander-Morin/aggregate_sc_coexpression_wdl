# syntax=docker/dockerfile:1

FROM rocker/r-ver:4.4.1 

# Install wget for conda installer and lipglpk40 is needed for Seurat
RUN apt-get update && apt-get install -y wget libglpk40 && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda (needed for anndata and reticulate)
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -afy

# Create conda environment (and accept TOS!) 
COPY docker/convert_h5ad_to_seurat/environment.yml .

RUN /opt/conda/bin/conda tos accept && \
    /opt/conda/bin/conda config --set channel_priority strict && \
    /opt/conda/bin/conda env create -f environment.yml

# Ensure anndata was installed correctly
RUN /opt/conda/envs/h5ad_env/bin/python -c 'import anndata'

# Install R packages (was getting errors with reticulate when using install.packages)
RUN R -q -e 'install.packages("remotes")' && \
    R -q -e 'remotes::install_github("rstudio/reticulate")' && \
    R -e "install.packages(c('Rtsne', 'sctransform', 'anndata', 'Seurat'), repos='https://cloud.r-project.org'))"

# Linking reticulate and python environment
ENV RETICULATE_CONDA="/opt/conda/bin/conda"
ENV RETICULATE_PYTHON="/opt/conda/envs/h5ad_env/bin/python"

# Get the conversion script and make it executable
COPY scripts/convert_h5ad_to_seurat.R /usr/local/bin/
RUN chmod +x /usr/local/bin/convert_h5ad_to_seurat.R

# Grab existing bioconductor docker image
FROM bioconductor/bioconductor_docker:RELEASE_3_21-r-4.5.1

# Was getting an error without installing KEGGREST first (why...)
RUN R -e 'BiocManager::install("KEGGREST")'
RUN R -e 'BiocManager::install("biomaRt")'

# Get the script that accesses biomart to download gene tables
COPY scripts/download_ensembl_pc_genes.R /usr/local/bin/

RUN chmod +x /usr/local/bin/download_ensembl_pc_genes.R
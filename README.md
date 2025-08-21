# aggregate_sc_coexpression_wdl

This project looks to containerize into a WDL framework parts of previous analysis that I performed here: https://github.com/PavlidisLab/TR_singlecell. This consisted of a series of R and bash scripts that required manual execution and environment setup, and I wanted to improve its reproducibility and portability.

Currently the worklow accepts a Cell x Gene URL and a given dataset identifer, and produces a gene-gene matrix of the aggregate coexpression scores, saved as a ```.tsv```. It was tested on a lightweight [dataset](https://cellxgene.cziscience.com/collections/0a77d4c0-d5d0-40f0-aa1a-5e1429bcbd7e).

### Workflow overview

In brief, the pipeline aims to:
1. Download scRNA-seq data (in ```.h5ad``` format) from the Cell x Gene resource
2. Convert the downloaded ```.h5ad``` file into a Seurat object
3. Preprocess this data by performing QC and filtering
4. Generate an aggregated coexpression matrix to assist in identifiyng gene-gene correlations that are reproducible across cell types

### Prereqs

- [Docker engine](https://docs.docker.com/engine/install/ubuntu/)
- Conda (I used Miniconda) 

``` wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh ```

This was developed and tested within a WSL2 environment (Ubuntu 22.04)

### Setup

1. Clone this repo
```
git clone https://github.com/Alexander-Morin/aggregate_sc_coexpression_wdl.git
cd aggregate_sc_coexpression_wdl
```
2. Create and actiate the Conda environment
```
conda env create -f environment.yml
conda activate wdl_env
```
3. Build the Docker images
```
# assumes you have started Docker with 'sudo service docker start'
bash scripts/build_images.sh
```

### Usage
1. Edit the input JSON (```inputs/main.json```) to specifiy your target data!
2. Run the workflow: ```miniwdl run workflows/main.wdl -i inputs/main.json```

The primary output will be 2 matrices each saved as ```.tsv``` files saved in the directory created by miniwdl:
1. The gene-gene matrix that scores the aggregated coexpression
2. A gene-gene matrix tracking how many times the gene pair WAS NOT measured together across cell types and thus did not produce a coexpression value to aggregate


# Convert a .h5ad file into a format usable by Seurat
# https://github.com/satijalab/seurat/issues/9072
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(anndata))
suppressPackageStartupMessages(library(Seurat))

use_python("/opt/conda/envs/h5ad_env/bin/python")
use_condaenv("h5ad_env")
py_config()

args <- commandArgs(trailingOnly = TRUE)

h5ad_path <- args[1]
outfile <- args[2]  # output file name (created in WDL)

message("Converting h5ad file to Seurat object")

h5ad <- read_h5ad(h5ad_path)

# Set up object using the raw counts
seurat_dat <- CreateSeuratObject(counts = t(as.matrix(h5ad$raw)), 
                                 meta.data = h5ad$obs,
                                 min.features = 0, # filtering downstream 
                                 min.cells = 0)

saveRDS(seurat_dat, outfile)

message(paste0("Successfully converted to Seurat object and saved to: ", outfile))
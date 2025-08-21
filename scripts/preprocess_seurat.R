# Load an .RDS of a Seurat object of raw counts and filters and normalizes it, 
# saving the cleaned data out as an .RDS
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Matrix))
source("preprocess_seurat_utils.R")

args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
seurat_raw_path <- args[2]
gene_table_path  <- args[3]
mt_table_path  <- args[4]
qc_plot1_path <- args[5]
qc_plot2_path <- args[6]
processed_file_path <- args[7]

message(paste0("Beginning to preprocess: ", id))

# Protein coding and mitochondrial gene tables
gene_table <- read.delim(gene_table_path , stringsAsFactors = FALSE)
mt_table <- read.delim(mt_table_path , stringsAsFactors = FALSE)

# Load raw seuarat data (over list generalizes to single or multi loading)
dat <- lapply(seurat_raw_path, readRDS)
dat <- reduce(dat, merge)
  
# Extract count matrix: default counts slot, but use data slot if counts empty
mat <- GetAssayData(dat, slot = "counts")
  
if (length(mat) == 0 || all(rowSums(mat) == 0)) {
  mat <- GetAssayData(dat, slot = "data")
}
  
stopifnot(typeof(mat@i) == "integer", "Expected that raw matrix has integer counts")
  
# Ready metadata
change_colnames <- c(Cell_type = "cell_type", Old_ID = "ID")
  
meta <- dat[[]] %>% 
  dplyr::rename(any_of(change_colnames)) %>% 
  rownames_to_column(var = "ID") %>% 
  add_count_info(mat = mat, mt_table = mt_table)
  
  
# QC plots
p1 <- all_hist(meta)
p2 <- qc_scatter(meta)
  
ggsave(p1, device = "png", dpi = 300, height = 12, width = 16, bg = "white",
       filename = qc_plot1_path)
  
ggsave(p2, device = "png", dpi = 300, height = 8, width = 8,
       filename = qc_plot2_path)
  
  
# Remove cells failing QC, keep only protein coding genes, and normalize
mat <- rm_low_qc_cells(mat, meta) %>%
  ensembl_to_symbol(ensembl_df = gene_table) %>% 
  get_pcoding_only(pcoding_df = gene_table) %>% 
  Seurat::NormalizeData(., normalization.method = "RC", scale.factor = 1e6, verbose = FALSE)
  
meta <- filter(meta, ID %in% colnames(mat))
mat <- mat[, meta$ID]
  
stopifnot(identical(colnames(mat), meta$ID), length(meta$ID) > 0)
  
message(paste("Count of cells:", ncol(mat),
                "Count unique cell types: ", n_distinct(meta$Cell_type)))
  
saveRDS(list(Mat = mat, Meta = meta), file = processed_file_path)

message(paste0(id, " successfully processed and saved to: ", processed_file_path))

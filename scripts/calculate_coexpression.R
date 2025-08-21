# Takes in a .RDS of a cleaned count matrix + metadata, performs coexpression
# on each cell type, and aggregates the result, saving out two .tsv files:
# the aggregated coexpression matrix, and a matrix tracking NAs/gene pairs
# that were not comeasured
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(aggtools))

args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
cleaned_data_path <- args[2]
gene_table_path  <- args[3]
agg_mat_path  <- args[4]
na_mat_path  <- args[5]

message(paste0("Beginning to calculate aggregate coexpression for: ", id))

gene_table <- read.delim(gene_table_path , stringsAsFactors = FALSE)

# Load processed data
dat <- readRDS(cleaned_data_path)
meta <- dat$Meta
mat <- dat$Mat

stopifnot(identical(colnames(mat), meta$ID))

# Run coexpression aggregate: returns a list with aggregate and NA tracking matrices
agg_l <- aggr_coexpr_single_dataset(
  mat = mat, 
  meta = meta,
  pc_df = gene_table,
  cor_method = "pearson",
  agg_method = "FZ"  # Averages Fisher's Z transformed correlation values
)

# Save output matrices as .tsv
fwrite_mat(agg_l$Agg_mat, agg_mat_path)
fwrite_mat(agg_l$NA_mat, na_mat_path)

message(paste0("Aggregate coexpression was successful for: ", id))

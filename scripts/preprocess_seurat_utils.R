# Functions for QC/preprocessing count matrices and metadata
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

# Gives a matrix with ENSEMBL IDs as rownames, return the matrix with the 
# corresponding gene symbols as rownames. Blank gene symbols are removed

ensembl_to_symbol <- function(mat, ensembl_df) {
  
  stopifnot(c("Gene_ID", "Symbol") %in% colnames(ensembl_df))
  
  ids <- intersect(pc$Gene_ID, rownames(mat))
  
  if (length(ids) == 0) stop("No common ENSEMBL IDs in rownames of matrix")
  
  common_genes <- data.frame(
    ID = ids,
    Symbol = ensembl_df$Symbol[match(ids, ensembl_df$Gene_ID)]) %>% 
    filter(Symbol != "")
    
  mat <- mat[common_genes$ID, ]
  rownames(mat) <- common_genes$Symbol
  
  return(mat)
}



# Assumes that mat is sparse gene x cell count matrix. Filters the matrix for 
# unique gene symbols in pcoding_df, and fills the missing genes as 0s. 

get_pcoding_only <- function(mat, pcoding_df) {
  
  stopifnot("Symbol" %in% colnames(pcoding_df))
  
  genes <- unique(pcoding_df$Symbol)
  common <- intersect(rownames(mat), genes)
  missing <- setdiff(genes, rownames(mat))
  
  if (length(common) == 0) stop("No common symbols in rownames of mat")
  
  pc_mat <- mat[common, ]
  pc_mat <- rbind(pc_mat, Matrix(0, nrow = length(missing), ncol = ncol(mat)))
  rownames(pc_mat) <- c(common, missing)
  pc_mat <- pc_mat[genes, ]
  
  return(pc_mat)
}



# Given a vector of genes that have either common gene symbols or ensembl
# ids, return a subst of gene_vec only containing the mitochondrial genes. 
# Assumes gene_vec has only mouse or human symbols/ensembl IDs.

get_mt_genes <- function(gene_vec, mt_table) {
  
  mt_genes <- gene_vec[gene_vec %in% c(mt_table$Gene_stable_ID, mt_table$Gene_name)]

  return(mt_genes)
}



# This adds columns to metadata: the number of total UMI counts for each 
# cell/column of mat, the number of non-zero expressing genes, and the RNA
# novelty/compexity, which is the ratio of the log10 gene counts to log10 umi 
# counts. It additionally adds ratio of mitochondrial if available. 
# https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html

add_count_info <- function(mat, meta, mt_table_path) {
  
  mt_genes <- get_mt_genes(rownames(mat), mt_table)
  
  meta <- meta %>% 
    mutate(
      UMI_counts = colSums(mat),
      Gene_counts = colSums(mat > 0),
      RNA_novelty = log10(Gene_counts) / log10(UMI_counts)
    )
  
  if (length(mt_genes) > 0) {
    mt_ratio <- colSums(mat[mt_genes, , drop = FALSE]) / meta$UMI_counts
    meta$MT_ratio = mt_ratio
  }
  
  # Remove cell x gene features of this type, if present
  meta <- meta[, !(colnames(meta) %in% c("nFeature_RNA", "nCount_RNA"))]
  
  return(meta)
}



# This subsets mat to remove cells that fail any of the filters laid out in: 
# https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html
# The RNA novelty filter is relaxed for Smart-seq runs as more reads can go to
# genes like mitochondrial https://pubmed.ncbi.nlm.nih.gov/33662621/

rm_low_qc_cells <- function(mat, 
                            meta,
                            min_counts = 500,
                            min_genes = 250,
                            min_novelty = NULL,
                            max_mt_ratio = 0.2) {
  
  keep_id <- meta %>%
    mutate(Is_smartseq = str_detect(str_to_lower(assay), "smart-seq")) %>%
    filter(
      UMI_counts >= min_counts,
      Gene_counts >= min_genes,
      ifelse(Is_smartseq, RNA_novelty > 0.5, RNA_novelty > 0.8)) %>%
    pull(ID)
  
  if ("MT_ratio" %in% colnames(meta)) {
    keep_id <- intersect(keep_id, filter(meta, MT_ratio < max_mt_ratio)$ID)
  }
  
  if (length(keep_id) == 0) stop("No remaining cells after filtering")
  
  return(mat[, keep_id])
}


# For QC plots for count matrices
# ------------------------------------------------------------------------------


qc_hist <- function(meta, xvar, xline, xlab, ylab, log10_xvar = TRUE) {
  
  if (log10_xvar) {
    meta[, xvar] <- log10(meta[, xvar])
    xline <- log10(xline)
  }
  
  ggplot(meta, aes(x = !!sym(xvar))) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = xline, colour = "red") +
    xlab(xlab) +
    ylab(ylab) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.margin = margin(10, 20, 10, 10))
}



all_hist <- function(meta) {
  
  pa <- qc_hist(meta, xvar = "UMI_counts", xline = 500, xlab = "Column (cell)-wise counts", ylab = "log10 UMI counts")
  pb <- qc_hist(meta, xvar = "Gene_counts", xline = 250, xlab = "Row (gene)-wise counts", ylab = "log10 Gene counts")
  pc <- qc_hist(meta, xvar = "RNA_novelty", xline = 0.8, xlab = "Column (cell)-wise counts", ylab = "log10 Gene counts / log10 UMI counts", log10_xvar = FALSE)
  
  p <- cowplot::plot_grid(pa, pb, pc, nrow = 2)
  
  if ("MT_ratio" %in% colnames(meta)) {
    pd <- qc_hist(meta, xvar = "MT_ratio", xline = 0.2, xlab = "Column (cell)-wise counts", ylab = "Ratio of mitochondrial counts", log10_xvar = FALSE)
    p <- cowplot::plot_grid(pa, pb, pc, pd, nrow = 2)
  }
  
  return(p)
}



qc_scatter <- function(meta) {
  
  ggplot(meta, aes(x = log10(Gene_counts), y = log10(UMI_counts))) +
    geom_point() +
    geom_vline(xintercept = log10(250), colour = "red") +
    geom_hline(yintercept = log10(500), colour = "red") +
    xlab("log10 Gene counts") +
    ylab("log10 UMI counts") +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20))
  
}

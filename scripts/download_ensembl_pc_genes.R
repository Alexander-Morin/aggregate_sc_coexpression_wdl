# Downloads human or mouse protein coding gene tables from ENSEMBL for specified
# ENSEMBL version
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(biomaRt))

args <- commandArgs(trailingOnly = TRUE)

species <- tolower(args[1])  # either human or mouse
version <- args[2]  # which ensembl version to use
outfile <- args[3]  # output file name (created in WDL)


if (species == "mouse") {
  
  symbol <- "mgi_symbol"
  species_data <- "mmusculus_gene_ensembl"
  chr_filter <- c(1:19, "MT", "X", "Y")
  
} else if (species == "human") {
  
  symbol = "hgnc_symbol"
  species_data = "hsapiens_gene_ensembl"
  chr_filter <- c(1:22, "MT", "X", "Y")
  
} else {
  
  stop("Species must be 'human' or 'mouse'")
  
}


download_ensembl_pcoding <- function(species_data,
                                     symbol,
                                     chr_filter,
                                     version) {

  # Which columns are acquired
  attributes <- c(
    "chromosome_name",
    "transcription_start_site",
    "transcript_start",
    "transcript_end",
    "strand",
    "ensembl_gene_id",
    symbol,
    "ucsc",
    "gene_biotype"
  )
  
  ens_mart <- useEnsembl(biomart = "ensembl",
                         dataset = species_data,
                         version = version)
  
  anno_table <- getBM(
    attributes = attributes,
    filters = "chromosome_name",
    values = chr_filter,
    mart = ens_mart,
    useCache = FALSE
  )
  
  # only protein coding gene type and order the table by chromosome then by TSS
  anno_table <- anno_table[anno_table$gene_biotype == "protein_coding", ]
  anno_table <- anno_table[order(match(anno_table$chromosome_name, chr_filter), anno_table$transcription_start_site), ]
  anno_table$gene_biotype <- NULL
  
  colnames(anno_table) <- c(
    "Chromosome",
    "Transcription_start_site",
    "Start",
    "End",
    "Strand",
    "Gene_ID",
    "Symbol",
    "Transcript_ID"
  )
  
 return(anno_table)
  
}


message(paste0("Downloading ENSEMBL protein coding table: ", 
               species, " version ", version))


anno_table <- download_ensembl_pcoding(species_data = species_data,
                                       symbol = symbol,
                                       chr_filter = chr_filter,
                                       version = version)

write.table(anno_table,
            quote = FALSE,
            row.names = FALSE,
            sep = "\t",
            file = outfile)

message("ENSEMBL protein coding table saved to ", outfile)

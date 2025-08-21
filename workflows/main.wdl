version 1.1

# Import tasks
import "./tasks/download_sc_data.wdl" as sc_data_downloader
import "./tasks/download_ensembl_pc_genes.wdl" as ensembl_gene_downloader
import "./tasks/convert_h5ad_to_seurat.wdl" as h5ad_to_seurat_converter
import "./tasks/preprocess_seurat.wdl" as seurat_preprocessor
import "./tasks/calculate_coexpression.wdl" as coexpression_calculator

# Main workflow: from downloading data to calculating aggregate coexpression
workflow main_wf {
  input {
    String dataset_id
    String dataset_url
    String species
    String ensembl_version
    File mt_gene_file
  }


  call sc_data_downloader.download_cxg_single_file as download_raw_data {
    input:
      dataset_id = dataset_id,
      dataset_url = dataset_url
  }


  call ensembl_gene_downloader.download_ensembl_pc_genes as get_ensembl_table {
    input:
      species = species,
      ensembl_version = ensembl_version
  }


  call h5ad_to_seurat_converter.convert_h5ad_to_seurat as h5ad_converter {
    input:
      dataset_id = dataset_id,
      h5ad_file = download_raw_data.h5ad_file
  }


  call seurat_preprocessor.preprocess_seurat as preprocessor {
    input:
      dataset_id = dataset_id,
      seurat_raw_file = h5ad_converter.seurat_raw_file,
      ensembl_table_file = get_ensembl_table.ensembl_file,
      mt_gene_file = mt_gene_file
  }


call coexpression_calculator.calculate_coexpression as coexpression {
    input:
      dataset_id = dataset_id
      cleaned_data_path = preprocessor.cleaned_mat_and_meta_file,
      ensembl_file_path = get_ensembl_table.ensembl_file
  }


}
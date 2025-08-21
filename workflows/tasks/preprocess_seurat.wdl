version 1.1

# Task: Produce a normalized and filtered count matrix and meta, saved as .RDS
# ---

task preprocess_seurat {
  input {
    String dataset_id
    File seurat_raw_file
    File ensembl_table_file
    File mt_gene_file
  }

  # Generate the output file names based on dataset_id
  String outfile_name = "~{dataset_id}_cleaned_mat_and_meta.RDS"
  String qc_plot1_name = "~{dataset_id}_scatter_plot.png"
  String qc_plot2_name = "~{dataset_id}_hist_plots.png"

  command {
    Rscript /usr/local/bin/preprocess_seurat.R \
      ${dataset_id} \
      ${seurat_raw_file} \
      ${ensembl_table_file} \
      ${mt_gene_file} \
      ${qc_plot1_name} \
      ${qc_plot2_name} \
      ${outfile_name}
  }

  output {
    File cleaned_mat_and_meta_file = outfile_name
    File qc_plot1 = qc_plot1_name
    File qc_plot2 = qc_plot2_name
  }

  runtime {
    docker: 'preprocess_seurat:1.0'  #  
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 20 SSD"
  }
}

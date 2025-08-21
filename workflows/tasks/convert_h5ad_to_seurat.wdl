version 1.1

# Task: Convert .h5ad file to a Seurat object saved as an .rds file
# ---

task convert_h5ad_to_seurat {
  input {
    String dataset_id
    File h5ad_file
  }

  # Generate the output file name based on dataset_id
  String outfile_name = "~{dataset_id}_Seurat_raw.RDS"

  command {
    Rscript /usr/local/bin/convert_h5ad_to_seurat.R \
      ${h5ad_file} \
      ${outfile_name}
  }

  output {
    File seurat_raw_file = outfile_name
  }

  runtime {
    docker: 'convert_h5ad:1.0'  #  
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 20 SSD"
  }
}

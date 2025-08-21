version 1.1

# Task: calculate aggregate coexpression for a dataset
# ---

task calculate_coexpression {
  input {
    String dataset_id
    File cleaned_data_file
    File ensembl_file
  }

  # Generate the output file name based on ID
  String agg_outfile_name = "aggregate_coexpression_mat_~{dataset_id}.tsv"
  String na_outfile_name = "na_tracking_mat_~{dataset_id}.tsv"

  command {
    Rscript /usr/local/bin/calculate_coexpression.R \
      ${dataset_id} \
      ${cleaned_data_file} \
      ${ensembl_file} \
      ${agg_outfile_name} \
      ${na_outfile_name}
  }

  output {
    File agg_coexpr_mat = agg_outfile_name
    File na_mat = na_outfile_name
  }

  runtime {
    docker: 'calculate_coexpression:1.0'  #  built on rocker/r-ver:4.4.1
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 20 SSD"
  }
}

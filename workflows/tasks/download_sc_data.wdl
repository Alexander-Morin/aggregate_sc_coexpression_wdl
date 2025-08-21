version 1.1

# Task: Download a single .h5ad file from Cell x Gene using wget
# ---

task download_cxg_single_file {
  input {
    String dataset_url
    String dataset_id
  }

  command <<<
    
    set -e -o pipefail
    
    # Download gets 3 tries, waiting 5 seconds between tries. File is saved
    # as <dataset_id>.h5ad, e.g. GSE12345.h5ad
    wget --quiet \
         --continue \
         --retry-connrefused \
         --waitretry=5 \
         --tries=3 \
         -O "~{dataset_id}.h5ad" \
         "~{dataset_url}"

    # Check for empty file
    if [ ! -s "~{dataset_id}.h5ad" ]; then
      echo "ERROR: Download resulted in an empty file."
      exit 1
    fi
  
  >>>

  output {
    File h5ad_file = "${dataset_id}.h5ad"
  }

  runtime {
    docker: 'download_sc_data:1.0'  # base image (ubuntu 22.04) with wget installed
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 20 SSD"
  }
}

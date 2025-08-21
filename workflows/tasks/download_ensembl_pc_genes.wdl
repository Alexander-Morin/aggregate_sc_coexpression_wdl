version 1.1

# Task: Download Ensembl Protein Coding Genes for mouse or human
# ---

task download_ensembl_pc_genes {
  input {
    String species
    String ensembl_version
  }

  # Generate the output file name based on species and version
  String outfile_name = "ensembl_protein_coding_~{species}_V~{ensembl_version}.tsv"

  command {
    Rscript /usr/local/bin/download_ensembl_pc_genes.R \
      ${species} \
      ${ensembl_version} \
      ${outfile_name}
  }

  output {
    File ensembl_file = outfile_name
  }

  runtime {
    docker: 'r-biomart:1.0'  #  bioconductor RELEASE_3_21-r-4.5.1 + biomart
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 20 SSD"
  }
}

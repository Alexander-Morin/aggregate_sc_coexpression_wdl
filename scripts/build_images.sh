#!/bin/bash

set -e

echo "Building Docker images..."

docker build -t download_sc_data:1.0 --network=host -f docker/download_sc_data/Dockerfile .
docker build -t r-biomart:1.0 --network=host -f docker/r-biomart/Dockerfile .
docker build -t convert_h5ad_to_seurat:1.0 --network=host -f docker/convert_h5ad_to_seurat/Dockerfile .
docker build -t preprocess_seurat:1.0 --network=host -f docker/preprocess_seurat/Dockerfile .
docker build -t calculate_coexpression:1.0 --network=host -f docker/calculate_coexpression/Dockerfile .

echo "Successfully built Docker images!"
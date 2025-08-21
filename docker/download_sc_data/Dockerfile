# syntax=docker/dockerfile:1

# Developed using WSL2 (Ubuntu 24.04.3)
FROM ubuntu:22.04

# Install wget and clean up the apt cache to keep the image small
RUN apt-get update && \
    apt-get install -y wget && \
    rm -rf /var/lib/apt/lists/*

# Set a working directory inside the container for data
WORKDIR /data
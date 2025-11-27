Docker First Pipeline
Overview

Docker First Pipeline is a bioinformatics workflow container for RNA-seq and ATAC-seq processing.

This repository includes:

run_rna.sh – RNA-seq processing pipeline

run_atac.sh – ATAC-seq processing pipeline with MACS2 peak calling

Dockerfile – builds a container with all required bioinformatics tools

test_data/ – fake FASTQ and genome files for testing CI workflows

GitHub Actions CI workflow to automatically test pipeline functionality

The pipelines are CI-friendly and can run end-to-end on fake/small datasets for testing.

Table of Contents

Installation

Running CI on GitHub

Running Locally

Test Data

Output Files

Notes

Installation

You need Docker installed to build and run the container.

Install Docker: https://docs.docker.com/get-docker/

Clone this repository:

git clone https://github.com/KhoiBio/docker_first_pipeline.git
cd docker_first_pipeline


Build the Docker container:

docker build -t docker_first_pipeline:latest .

Running CI on GitHub

This repository includes a GitHub Actions workflow that automatically tests the pipelines using fake small data.

Workflow location: .github/workflows/ci.yml

Triggered on push or pull requests to main branch

The CI workflow performs:

Docker build

Runs RNA-seq pipeline on test_data/fake_rna_R1.fastq.gz & R2

Runs ATAC-seq pipeline on test_data/fake_atac_R1.fastq.gz & R2

Checks for key output files

Expected CI outputs:

fake_rna_final.bam for RNA-seq

condition1_merged.sorted.bam for ATAC-seq

You can see CI results in the Actions tab of your GitHub repository.

Running Locally

To test the pipelines on your local machine:

RNA-seq pipeline
docker run --rm \
  -v $PWD/test_data:/data \
  docker_first_pipeline:latest \
  bash /pipeline/workflow/run_rna.sh \
  /data/fake_rna_R1.fastq.gz \
  /data/fake_rna_R2.fastq.gz \
  /data/fake_genome.fa \
  /data/output

ATAC-seq pipeline
docker run --rm \
  -v $PWD/test_data:/data \
  docker_first_pipeline:latest \
  bash /pipeline/workflow/run_atac.sh \
  /data/fake_atac_R1.fastq.gz \
  /data/fake_atac_R2.fastq.gz \
  /data/fake_genome.fa \
  condition1 \
  /data/output

Notes:

The /data/output directory will contain all output files.

You can replace fake_* files with real FASTQ and reference genome for actual analysis.

Test Data

The repository contains fake small files for CI testing:

test_data/fake_rna_R1.fastq.gz & R2 – fake RNA FASTQ

test_data/fake_atac_R1.fastq.gz & R2 – fake ATAC FASTQ

test_data/fake_genome.fa – fake reference genome

These allow full pipeline execution in GitHub Actions or locally without large data.

Output Files

After running the pipelines, the key output files are:

Pipeline	Key Output Files
RNA-seq	final/<sample>_final.bam
ATAC-seq	merge/<condition>_merged.sorted.bam
macs2/all/<condition>_peaks.narrowPeak
Notes

The pipelines are designed to be CI-friendly and work with tiny fake datasets.

When using real data, ensure enough memory and CPU cores are available.

Docker ensures all dependencies and tools are contained, providing reproducible results.

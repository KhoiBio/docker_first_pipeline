
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Use bash for RUN commands
SHELL ["/bin/bash", "-c"]

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget curl unzip git build-essential python3 python3-pip \
    openjdk-11-jre-headless \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    libncurses5-dev libncursesw5-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and install Python packages
RUN pip3 install --upgrade pip && \
    pip3 install pandas numpy matplotlib seaborn MACS2

# Install bioinformatics tools via apt
RUN apt-get update && apt-get install -y \
    fastqc \
    samtools \
    bedtools \
    bowtie2 \
    hisat2 \
    subread \
    trim-galore \
    bwa \
    bamtools \
    && rm -rf /var/lib/apt/lists/*

# Install STAR aligner (RNA-seq)
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz && \
    tar -xvzf 2.7.10b.tar.gz && \
    cd STAR-2.7.10b/source && make STAR && cp STAR /usr/local/bin && \
    cd / && rm -rf STAR-2.7.10b 2.7.10b.tar.gz

# Install Picard (Java-based)
RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar -O /usr/local/bin/picard.jar

# Set working directory
WORKDIR /pipeline/workflow

# Copy pipeline scripts into container
COPY workflow/run_rna.sh workflow/run_atac.sh ./

# Make scripts executable
RUN chmod +x run_rna.sh run_atac.sh

# Set entrypoint
ENTRYPOINT ["/bin/bash","-c"]

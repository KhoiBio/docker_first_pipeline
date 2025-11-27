FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget curl unzip git build-essential python3 python3-pip \
    openjdk-11-jre-headless \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    libncurses5-dev libncursesw5-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install bioinformatics tools
RUN pip3 install --upgrade pip && \
    pip3 install pandas numpy matplotlib seaborn

# Tools for RNA/ATAC-seq
RUN apt-get update && apt-get install -y \
    fastqc \
    samtools \
    bedtools \
    bowtie2 \
    hisat2 \
    subread \
    trim-galore \
    && rm -rf /var/lib/apt/lists/*

# Install STAR aligner (RNA-seq)
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz && \
    tar -xvzf 2.7.10b.tar.gz && \
    cd STAR-2.7.10b/source && make STAR && cp STAR /usr/local/bin

# Install MACS2 (ATAC-seq)
RUN pip3 install MACS2

# Set working directory
WORKDIR /pipeline

ENTRYPOINT ["/bin/bash"]


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
    && rm -rf /var/lib/apt/lists/*

# Install STAR aligner (RNA-seq)
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz && \
    tar -xvzf 2.7.10b.tar.gz && \
    cd STAR-2.7.10b/source && make STAR && cp STAR /usr/local/bin && \
    cd / && rm -rf STAR-2.7.10b 2.7.10b.tar.gz

# Install BWA (ATAC-seq)
RUN apt-get update && apt-get install -y bwa && rm -rf /var/lib/apt/lists/*

# Install Picard (Java-based)
RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar -O /usr/local/bin/picard.jar

# Install BAMTools
RUN apt-get update && apt-get install -y bamtools && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /pipeline/workflow

# Copy pipeline scripts into container
COPY run_rna.sh run_atac.sh run_atac_macs.sh ./

# Make scripts executable
RUN chmod +x run_rna.sh run_atac.sh run_atac_macs.sh

# Set entrypoint
ENTRYPOINT ["/bin/bash"]

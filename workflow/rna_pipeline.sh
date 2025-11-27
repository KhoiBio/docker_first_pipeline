#!/bin/bash

FASTQ=$1  # e.g., sample_R1.fastq.gz
GENOME_INDEX=$2  # STAR genome index directory
GTF=$3  # annotation file (GTF)
OUTDIR="rna_output"
SAMPLE=$(basename "$FASTQ" | cut -d_ -f1)

mkdir -p $OUTDIR

# Step 1: Trim
trim_galore $FASTQ 

# Step 2: Align with STAR
STAR --genomeDir $GENOME_INDEX \
     --readFilesIn ${SAMPLE}_RNA_trimmed.fq \
     --readFilesCommand zcat \
     --outFileNamePrefix $OUTDIR/${SAMPLE}_ \
     --outSAMtype BAM SortedByCoordinate

# Step 3: Quantification
featureCounts -a $GTF -o $OUTDIR/${SAMPLE}_counts.txt $OUTDIR/${SAMPLE}_Aligned.sortedByCoord.out.bam

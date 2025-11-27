#!/usr/bin/bash
FASTQ=$1  # e.g., sample_R1.fastq.gz
GENOME_INDEX=$3  # Bowtie2 index basename
OUTDIR="atac_output"
SAMPLE=$(basename "$FASTQ" | cut -d_ -f1)

mkdir -p $OUTDIR

# Step 1: Trim
trim_galore --paired $FASTQ ${FASTQ/R1/R2}

# Step 2: Align
mv ${SAMPLE}_chr22_enriched_R2_val_2.fq.gz ${SAMPLE}_R2_val_2.fq.gz
mv ${SAMPLE}_chr22_enriched_R1_val_1.fq.gz ${SAMPLE}_R1_val_1.fq.gz


bowtie2 -x $GENOME_INDEX \
  -1 ${SAMPLE}_R1_val_1.fq.gz \
  -2 ${SAMPLE}_R2_val_2.fq.gz \
| samtools view -bS - \
| samtools sort -o $OUTDIR/${SAMPLE}_sorted.bam

samtools index $OUTDIR/${SAMPLE}_sorted.bam

# Step 4: Peak calling
macs2 callpeak -t $OUTDIR/${SAMPLE}_sorted.bam -f BAMPE -n $SAMPLE -g hs --outdir $OUTDIR

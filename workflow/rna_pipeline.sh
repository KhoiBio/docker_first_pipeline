
#!/bin/bash
set -e  # Exit on error

# Positional arguments
R1=$1        # First FASTQ file
R2=$2        # Second FASTQ file
GENOME=$3    # Genome index (e.g., /data/reference/hg38.fa)

# Optional: output directory and blacklist as arguments or defaults
OUTPUT_DIR=${4:-"/data/output"}
BLACKLIST=${5:-"/data/reference/GRCh38_unified_blacklist.bed"}

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate atac_env

# Derive sample name from R1 filename
SAMPLE=$(basename "$R1" | sed -E 's/_r1\.fastq\.gz//')

echo "Running ATAC pipeline for sample: $SAMPLE"
mkdir -p ${OUTPUT_DIR}/{trimmed,aligned,filtered,final}

# Trim adapters
trim_galore --paired --nextera --cores 12 -o ${OUTPUT_DIR}/trimmed $R1 $R2

TRIM_R1="${OUTPUT_DIR}/trimmed/${SAMPLE}_r1_val_1.fq.gz"
TRIM_R2="${OUTPUT_DIR}/trimmed/${SAMPLE}_r2_val_2.fq.gz"

# Align
bwa mem -t 12 $GENOME $TRIM_R1 $TRIM_R2 | samtools view -bS - > ${OUTPUT_DIR}/aligned/${SAMPLE}.bam

# Sort and index
samtools sort -@ 8 -o ${OUTPUT_DIR}/aligned/${SAMPLE}_sorted.bam ${OUTPUT_DIR}/aligned/${SAMPLE}.bam
samtools index ${OUTPUT_DIR}/aligned/${SAMPLE}_sorted.bam

# Remove duplicates
picard MarkDuplicates \
    I=${OUTPUT_DIR}/aligned/${SAMPLE}_sorted.bam \
    O=${OUTPUT_DIR}/filtered/${SAMPLE}_dedup.bam \
    M=${OUTPUT_DIR}/filtered/${SAMPLE}_dedup_metrics.txt \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT

samtools index ${OUTPUT_DIR}/filtered/${SAMPLE}_dedup.bam

# Remove mitochondrial reads
samtools view -h ${OUTPUT_DIR}/filtered/${SAMPLE}_dedup.bam \
    | grep -E '^@|^chr([1-9]|1[0-9]|2[0-2]|X|Y)\b' \
    | samtools view -b -o ${OUTPUT_DIR}/filtered/${SAMPLE}_nomt.bam

# Remove blacklisted regions
bedtools intersect -v -abam ${OUTPUT_DIR}/filtered/${SAMPLE}_nomt.bam -b $BLACKLIST > ${OUTPUT_DIR}/filtered/${SAMPLE}_noblacklisted.bam

# Filter flags
samtools view -b -F 0x904 ${OUTPUT_DIR}/filtered/${SAMPLE}_noblacklisted.bam > ${OUTPUT_DIR}/filtered/${SAMPLE}_filtered.bam

# BAMTools filters
bamtools filter -in ${OUTPUT_DIR}/filtered/${SAMPLE}_filtered.bam -out ${OUTPUT_DIR}/filtered/${SAMPLE}_bamtools.bam \
    -tag "NM:<=4" \
    -isProperPair true \
    -insertSize "<=2000" \
    -isMapped true \
    -isPrimaryAlignment true \
    -isDuplicate false \
    -isPaired true \
    -isMateMapped true

# Final filter for proper pair and same chromosome
samtools view -h ${OUTPUT_DIR}/filtered/${SAMPLE}_bamtools.bam \
    | awk 'BEGIN {OFS="\t"} /^@/ {print; next} $7 == "=" && and($2, 0x2) {print}' \
    | samtools view -b -o ${OUTPUT_DIR}/final/${SAMPLE}_final.bam

samtools index ${OUTPUT_DIR}/final/${SAMPLE}_final.bam

echo "Pipeline completed for $SAMPLE"

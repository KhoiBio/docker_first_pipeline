
#!/bin/bash
set -e

# Positional arguments
R1=$1
R2=$2
GENOME=$3
CONDITION=$4
OUTPUT_DIR=${5:-"/data/output"}
BLACKLIST=${6:-"/data/reference/GRCh38_unified_blacklist.bed"}
GENOME_SIZE=${7:-2913022398}

eval "$(conda shell.bash hook)"
conda activate atac_env

SAMPLE=$(basename "$R1" | sed -E 's/_r1\.fastq\.gz//')

echo "Running ATAC + MACS2 pipeline for sample: $SAMPLE (Condition: $CONDITION)"
mkdir -p ${OUTPUT_DIR}/{trimmed,aligned,filtered,final,merge,macs2/all}

# Step 1: Trim adapters
trim_galore --paired --nextera --cores 12 -o ${OUTPUT_DIR}/trimmed $R1 $R2
TRIM_R1="${OUTPUT_DIR}/trimmed/${SAMPLE}_r1_val_1.fq.gz"
TRIM_R2="${OUTPUT_DIR}/trimmed/${SAMPLE}_r2_val_2.fq.gz"

# Step 2: Align
bwa mem -t 12 $GENOME $TRIM_R1 $TRIM_R2 | samtools view -bS - > ${OUTPUT_DIR}/aligned/${SAMPLE}.bam

# Step 3: Sort and index
samtools sort -@ 8 -o ${OUTPUT_DIR}/aligned/${SAMPLE}_sorted.bam ${OUTPUT_DIR}/aligned/${SAMPLE}.bam
samtools index ${OUTPUT_DIR}/aligned/${SAMPLE}_sorted.bam

# Step 4: Remove duplicates
picard MarkDuplicates \
    I=${OUTPUT_DIR}/aligned/${SAMPLE}_sorted.bam \
    O=${OUTPUT_DIR}/filtered/${SAMPLE}_dedup.bam \
    M=${OUTPUT_DIR}/filtered/${SAMPLE}_dedup_metrics.txt \
    REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
samtools index ${OUTPUT_DIR}/filtered/${SAMPLE}_dedup.bam

# Step 5: Remove mitochondrial reads
samtools view -b ${OUTPUT_DIR}/filtered/${SAMPLE}_dedup.bam \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
    > ${OUTPUT_DIR}/filtered/${SAMPLE}_nomt.bam

# Step 6: Remove blacklist
bedtools intersect -v -abam ${OUTPUT_DIR}/filtered/${SAMPLE}_nomt.bam -b $BLACKLIST > ${OUTPUT_DIR}/filtered/${SAMPLE}_noblacklisted.bam

# Step 7: Filter flags
samtools view -b -F 0x904 ${OUTPUT_DIR}/filtered/${SAMPLE}_noblacklisted.bam > ${OUTPUT_DIR}/filtered/${SAMPLE}_filtered.bam

# Step 8: BAMTools filter
bamtools filter -in ${OUTPUT_DIR}/filtered/${SAMPLE}_filtered.bam -out ${OUTPUT_DIR}/filtered/${SAMPLE}_bamtools.bam \
    -tag "NM:<=4" -isProperPair true -insertSize "<=2000" -isMapped true -isPrimaryAlignment true -isDuplicate false -isPaired true -isMateMapped true

# Step 9: Final BAM
samtools view -h ${OUTPUT_DIR}/filtered/${SAMPLE}_bamtools.bam \
    | awk 'BEGIN {OFS="\t"} /^@/ {print; next} $7 == "=" && and($2, 0x2) {print}' \
    | samtools view -b -o ${OUTPUT_DIR}/final/${SAMPLE}_final.bam
samtools index ${OUTPUT_DIR}/final/${SAMPLE}_final.bam

# Step 10: Merge replicates by condition
samtools merge ${OUTPUT_DIR}/merge/${CONDITION}_merged.bam ${OUTPUT_DIR}/final/${CONDITION}_*_final.bam
samtools sort -o ${OUTPUT_DIR}/merge/${CONDITION}_merged.sorted.bam ${OUTPUT_DIR}/merge/${CONDITION}_merged.bam
samtools index ${OUTPUT_DIR}/merge/${CONDITION}_merged.sorted.bam
rm ${OUTPUT_DIR}/merge/${CONDITION}_merged.bam

# Step 11: Convert to BEDPE
samtools sort -n -o ${OUTPUT_DIR}/merge/${CONDITION}_merged.name_sorted.bam ${OUTPUT_DIR}/merge/${CONDITION}_merged.sorted.bam
bedtools bamtobed -bedpe -i ${OUTPUT_DIR}/merge/${CONDITION}_merged.name_sorted.bam > ${OUTPUT_DIR}/merge/${CONDITION}_merged.bedpe

# Step 12: Shift Tn5 sites
awk -F $'\t' 'BEGIN {OFS = FS}{ if ($9 == "+") {$2 = $2 + 4; $6 = $6 - 5} else if ($9 == "-") {$3 = $3 - 5; $5 = $5 + 4} print $0}' \
    ${OUTPUT_DIR}/merge/${CONDITION}_merged.bedpe > ${OUTPUT_DIR}/merge/${CONDITION}_merged.shifted.bedpe

# Step 13: MACS2 peak calling
macs2 callpeak -t ${OUTPUT_DIR}/merge/${CONDITION}_merged.shifted.bedpe \
    -f BEDPE -g ${GENOME_SIZE} -n ${CONDITION} \
    --outdir ${OUTPUT_DIR}/macs2/all --keep-dup all --call-summits


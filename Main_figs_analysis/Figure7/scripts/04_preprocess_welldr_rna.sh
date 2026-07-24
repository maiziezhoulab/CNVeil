#!/bin/bash
set -euo pipefail

# Adapted from: https://github.com/navinlabcode/wellDR-seq

RAW_FASTQ_DIR="example_data/figure7/transcriptomics/raw_fastq"
OUTPUT_DIR="results/figure7/transcriptomics"
WELLDR_CONVERTER="external/wellDR-seq/rna_fq_preprocess/scDR-fq-convert-v4.pl"
BARCODE_FILE="external/wellDR-seq/rna_fq_preprocess/com_index_wafer.txt"
TRIMMOMATIC_JAR="external/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTER_FILE="external/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"
STAR_INDEX="external/reference/hg19_star"
THREADS=20
SAMPLES=(SRR28357500 SRR28357501 SRR28357502)

mkdir -p "$OUTPUT_DIR"/{demux_logs,trimmed,star_out,demultiplexed_fastq}

for sample in "${SAMPLES[@]}"; do
  perl "$WELLDR_CONVERTER" \
    "$RAW_FASTQ_DIR/${sample}_1.fastq.gz" \
    "$RAW_FASTQ_DIR/${sample}_2.fastq.gz" \
    "$BARCODE_FILE" 20 \
    > "$OUTPUT_DIR/demux_logs/${sample}.log" 2>&1
done

CELL_FASTQ_DIR="$OUTPUT_DIR/demultiplexed_fastq"

for fastq in "$CELL_FASTQ_DIR"/*.fq.gz; do
  name=$(basename "$fastq" .fq.gz)
  java -jar "$TRIMMOMATIC_JAR" SE \
    -threads 1 \
    "$fastq" \
    "$OUTPUT_DIR/trimmed/${name}_trimmed.fq.gz" \
    "ILLUMINACLIP:${ADAPTER_FILE}:2:30:10" \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
done

MANIFEST="$OUTPUT_DIR/star_out/manifest.tsv"
: > "$MANIFEST"
for fastq in "$OUTPUT_DIR"/trimmed/*_trimmed.fq.gz; do
  cell=$(basename "$fastq" _trimmed.fq.gz)
  printf "%s\t-\t%s\n" "$fastq" "$cell" >> "$MANIFEST"
done

STAR \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUTPUT_DIR/star_out/wellDR_" \
  --soloType SmartSeq \
  --readFilesManifest "$MANIFEST" \
  --soloUMIdedup Exact NoDedup \
  --soloStrand Unstranded \
  --soloFeatures Gene GeneFull \
  --outSAMtype BAM SortedByCoordinate

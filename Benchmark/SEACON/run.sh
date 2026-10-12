#!/usr/bin/env bash
set -euo pipefail
# Activate SEACON before running this script.
seacon_output="/path/to/SEACON/TN1"
ref="/path/to/hg19.fa"
dataset_bams="/path/to/individual_cell_bams"
phased_vcf="/path/to/concatenated.vcf.gz"
matched_normal_bam="/path/to/hg19_matched_normal.bam"
map_file="/path/to/hg19_mappability.bigWig"
mkdir -p "$seacon_output"
seacon prep_readcount -o "$seacon_output" -i "$dataset_bams" -r "$ref" \
    -b 500000 --map-file "$map_file" -P 12 --normal-bam "$matched_normal_bam"
seacon prep_baf -o "$seacon_output" -i "$dataset_bams" --vcf "$phased_vcf" \
    --block-size 50 --normal-bam "$matched_normal_bam"
seacon call -o "$seacon_output" --tolerance 0.2 --upper-filter 5 --max-wgd 2

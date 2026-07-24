#!/bin/bash
set -euo pipefail

HAPLOTYPE_INPUT="example_data/figure7/copy_number/input_medicc2.tsv"
TOTAL_CN_INPUT="example_data/figure7/copy_number/subclone_totalcn.tsv"
HAPLOTYPE_OUT="results/figure7/medicc2_haplotype"
TOTAL_CN_OUT="results/figure7/medicc2_total_cn"
THREADS=10
BOOTSTRAPS=100

mkdir -p "$HAPLOTYPE_OUT" "$TOTAL_CN_OUT"

medicc2 "$HAPLOTYPE_INPUT" "$HAPLOTYPE_OUT" \
  --plot auto \
  --events \
  --bootstrap-method chr-wise \
  --bootstrap-nr "$BOOTSTRAPS" \
  -j "$THREADS"

medicc2 "$TOTAL_CN_INPUT" "$TOTAL_CN_OUT" \
  --total-copy-numbers \
  --input-allele-columns total_cn \
  --plot auto \
  --events \
  --bootstrap-method chr-wise \
  --bootstrap-nr "$BOOTSTRAPS" \
  -j "$THREADS"

#!/usr/bin/env bash
set -euo pipefail
source ~/.bashrc
CHISEL_ENV="chisel"
SPRINTER_ENV="sprinter"
REF="/path/to/hg19/genome.fa"
TUM="/path/to/all_tumor_cells.bam"
OUTPUT_DIR="/path/to/SPRINTER/TN1"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
conda activate "$CHISEL_ENV"
chisel_rdr -t "$TUM" -r "$REF" -b 50000 -m 100000 -j 8
gzip -c rdr/rdr.tsv > sprinter.input.tsv.gz
# Original deletion was optional housekeeping; retain the intermediate files.
conda deactivate
conda activate "$SPRINTER_ENV"
sprinter ./sprinter.input.tsv.gz --refgenome "$REF" --minreads 100000 \
    --rtreads 200 --cnreads 1000 --minnumcells 15 -j 8

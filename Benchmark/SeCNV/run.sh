#!/usr/bin/env bash
set -euo pipefail
SCRIPTS_DIR="/path/to/SeCNV/Scripts"
INPUT_DIR="/path/to/hg19_individual_cell_bams"
OUTPUT_DIR="/path/to/SeCNV/TN1"
REF="/path/to/hg19.fa"
cd "$SCRIPTS_DIR"
python SeCNV.py "$INPUT_DIR" "$OUTPUT_DIR" "$REF" --bin_size=50000

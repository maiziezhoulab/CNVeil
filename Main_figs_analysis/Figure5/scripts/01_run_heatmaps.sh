#!/bin/bash

PYTHON=python2.7
SCRIPT="scripts/01_generate_heatmaps.py"

# List of datasets to process
DATASETS=(
    S0_section_E
    TN1
    TN3
    TN4
    #...
)

for DATASET in "${DATASETS[@]}"; do
    echo "Processing ${DATASET}..."

    ${PYTHON} ${SCRIPT} \
        --folder example_data/figure5/${DATASET} \
        --barcode_file example_data/barcodes/${DATASET}_barcodes.txt \
        --cnveil_tsv example_data/figure5/${DATASET}/CNVeil.tsv \
        --pattern "*.tsv" \
        --format png \
        --gridsize 12,6
done
 

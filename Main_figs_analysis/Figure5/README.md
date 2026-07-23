# Figure 5: Allele-specific copy number comparison

This directory contains the scripts used to generate Figure 5 of the CNVeil manuscript.

## Panels

- **Panels a–d:** Genome-wide allele-specific copy number heatmaps for 10X S0 Section E, TN1, TN3, and TN4.
- **Panels e–j:** B-allele frequency distributions and allele-specific copy number profiles for selected genomic regions.

## Input

Required example files are stored under:

```text
example_data/figure5/
```

Each dataset folder contains allele-specific copy number results from:

- CNVeil
- CHISEL
- Alleloscope
- CNRein

A barcode mapping file is also required to match cell identifiers across methods.

## Workflow

### 1. Generate and reorder heatmaps

```bash
bash scripts/01_run_heatmaps.sh
```

This step:

- converts cell barcodes to consistent SRR identifiers;
- determines the cell order from CNVeil;
- applies the same order to all methods;
- generates one allele-specific copy number heatmap per method.

### 2. Generate selected-region panels

```bash
python scripts/02_plot_selected_regions.py
```

This step generates the BAF histograms and regional allele-specific copy number comparisons used in panels e–j.

## Output

The workflow generates:

- Figure 5a–d genome-wide heatmaps
- Figure 5e–j selected-region comparisons

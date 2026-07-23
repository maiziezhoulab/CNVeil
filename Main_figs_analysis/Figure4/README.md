# Ploidy Summary Analysis

This workflow summarizes genome-wide ploidy estimates across ACT and 10X datasets and generates the benchmark visualizations used in multiple figures of the CNVeil manuscript.

## Figures

This workflow generates the following manuscript panels:

- **Figure 4a** — Ploidy estimates compared with FACS measurements across ACT samples.
- **Figure S3** — Heatmap of ploidy estimates across ACT samples and methods.
- **Figure 3c** — Overall weighted MSE comparison on the T10 dataset.

## Input

Required files are available under `example_data/`.

Main inputs include:

- Ploidy summary table
- ACT ploidy reference (FACS)
- T10 benchmark metrics
- Overall weighted MSE summary

## Script

Run

```bash
python scripts/benchmark_ploidy_summary.py
```

## Output

This workflow produces:

- Figure 4a
- Figure S3
- Figure 3c

along with intermediate summary tables used for plotting.

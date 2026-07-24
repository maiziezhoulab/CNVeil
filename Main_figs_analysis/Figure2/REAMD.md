# Figure 2: Simulation benchmark

This directory contains the workflow used to evaluate total copy-number inference on simulated single-cell DNA sequencing datasets.

## Input

Example inputs are provided under:

```text
example_data/figure2/
```

Main inputs include:

- simulated ground-truth CN profiles;
- predicted cell-by-bin CN matrices from each tool;
- benchmark breakpoint annotations;
- genome FASTA index;
- plotting configuration table.

## Workflow

### 1. Evaluate breakpoint and CN-state accuracy

```bash
python scripts/01_evaluate_breakpoints.py \
  --comp_file example_data/figure2/predictions/CNVeil_ploidy5.csv \
  --bench_dir example_data/figure2/benchmark/ploidy5 \
  --output_dir results/figure2/CNVeil/ploidy5 \
  --max_dist_thresh 2000000 \
  --fai_path example_data/reference/genome.fa.fai
```

This step compares predicted and ground-truth breakpoints using:

- segment-only matching;
- CNV-aware matching requiring both breakpoint position and CN transition agreement.

### 2. Generate benchmark violin plots

Prepare a configuration file containing:

```text
tool,ploidy,eval_dir
```

Then run:

```bash
python scripts/02_plot_benchmark_violins.py \
  --input_path example_data/figure2/eval_config.csv \
  --output_dir results/figure2 \
  --max_dist_thresh 2000000
```

This step generates the recall, precision, and F1 distributions shown in Figure 2a.

### 3. Generate representative heatmaps

```bash
python scripts/03_plot_representative_heatmaps.py
```

This step plots the ground-truth and inferred integer CN profiles for the representative ploidy-5 simulation shown in Figure 2b.

## Output

The workflow reproduces:

- **Figure 2a:** segmentation and integer CN benchmark distributions;
- **Figure 2b:** representative ground-truth and inferred CN heatmaps.

# Figure 6: Bulk and pseudo-bulk copy number concordance

This directory contains the workflow used to compare CNVeil pseudo-bulk copy number profiles with matched bulk WES-derived ASCAT profiles.

The same workflow is also used for Figure S7 across the remaining ACT samples.

## Prerequisites

Bulk WES copy number profiles are generated using ASCAT.

Please install ASCAT and prepare the required reference files following the official documentation:

- :contentReference[oaicite:0]{index=0}
- :contentReference[oaicite:1]{index=1}

The detailed instructions for preparing targeted sequencing reference files, running allele counting, and generating ASCAT copy number segments are maintained by the ASCAT developers and are not duplicated in this repository. :contentReference[oaicite:2]{index=2}


## Input

Required example files are stored under:

```text
example_data/figure6/
```

Main inputs include:

- Matched tumor and normal WES BAM files
- ASCAT worksheet and targeted SNP reference files
- ASCAT segment files
- CNVeil pseudo-bulk haplotype copy number profiles
- ACT sample information table

## Workflow

### 1. Run ASCAT

```bash
Rscript scripts/01_run_ascat.R \
  <sample_id> \
  <normal_id> \
  <tumor_id> \
  <tumor_bam> \
  <normal_bam> \
  <gender> \
  <worksheet>
```

This step generates the bulk WES-derived ASCAT copy number segments.

### 2. Compare ASCAT and CNVeil pseudo-bulk profiles

```bash
python scripts/02_plot_bulk_pseudobulk_concordance.py
```

This step:

- maps ASCAT segments to the CNVeil pseudo-bulk bin grid;
- aligns the two bulk haplotypes with the CNVeil haplotypes;
- compares total and haplotype-specific copy number profiles;
- computes per-chromosome Pearson correlations;
- generates the combined 05–06–07 figure.

## Output

The workflow generates:

- **Figure 6** for TN6
- **Figure S7** for the remaining ACT samples
- Total copy number comparison tracks
- Allele1 and allele2 comparison tracks
- Per-chromosome correlation bar plots

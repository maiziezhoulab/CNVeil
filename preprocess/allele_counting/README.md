# VarTrix Allele Matrix Pipeline

This repository provides a reproducible pipeline to extract and process **allele-specific count matrices** from single-cell sequencing data using VarTrix.

The workflow consists of two main steps:

1. Extract per-chromosome REF/ALT matrices from BAM + VCF
2. Combine and preprocess matrices into genome-wide datasets

This pipeline is designed for downstream analyses such as:

* allele-specific copy number inference
* clonal structure reconstruction
* tumor heterogeneity analysis

---

## Overview

Given:

* a **10x Genomics BAM file**
* **phased VCF files (per chromosome)**
* a **cell barcode list**

the pipeline produces:

* genome-wide REF allele matrix
* genome-wide ALT allele matrix
* filtered SNP annotation table

---

## Repository Structure

```text
.
├── README.md
├── scripts/
│   ├── run_vartrix.sh
│   └── combine_vartrix_matrices.R
├── example/       
└── output/         
```

---

## Requirements

### Core tools

* VarTrix
* samtools
* tabix

### R packages

```r
install.packages("Matrix")
install.packages("optparse")
```

---

## Input Requirements

### 1. BAM file

10x Genomics BAM file:

```text
possorted_bam.bam
possorted_bam.bam.bai
```

---

### 2. Phased VCF files

Expected naming:

```text
<sample>_phased_chr1.vcf.gz
<sample>_phased_chr2.vcf.gz
...
```

Each must be:

* bgzipped
* indexed (`.tbi`)

---

### 3. Barcode file

```text
AAACCTGAGAAACCAT-1
AAACCTGAGAAACGAG-1
...
```

---

## Step 1: Run VarTrix

Extract REF/ALT matrices per chromosome.

```bash
bash scripts/run_vartrix.sh \
  --patient A \
  --bam /path/to/possorted_bam.bam \
  --ref /path/to/genome.fa \
  --barcodes /path/to/barcodes.txt \
  --vcf-dir /path/to/phased_vcfs \
  --out-dir ./output \
  --vartrix /path/to/vartrix_linux
```

### Output

```text
output/
├── chr1_matrix_ref
├── chr1_matrix_alt
├── chr2_matrix_ref
├── chr2_matrix_alt
...
```

---

## Step 2: Combine Matrices

Merge per-chromosome matrices into genome-wide matrices.

```bash
Rscript scripts/combine_vartrix_matrices.R \
  --dir ./output \
  --vcf-dir ./vcf_sub \
  --barcodes ./barcodes.tsv \
  --sample A \
  --out-dir ./output \
  --plot-qc
```

---

## Final Output

```text
output/
├── alt_all.rds
├── ref_all.rds
├── var_all.rds
├── A_vartrix_qc.pdf
```

### Description

* `alt_all.rds` → ALT allele count matrix
* `ref_all.rds` → REF allele count matrix
* `var_all.rds` → SNP metadata
* `*_qc.pdf` → QC plots (coverage + VAF)

---

## Filtering Strategy

Variants are filtered by:

```r
rowSums(alt) > 0
```

This removes SNPs with no ALT support across all cells.

---

## Workflow Summary

```text
BAM + VCF + barcodes
        ↓
   (VarTrix)
        ↓
chr-wise REF/ALT matrices
        ↓
 (R postprocessing)
        ↓
Genome-wide matrices (RDS)
```

---

## Example

```bash
bash scripts/run_vartrix.sh \
  --patient demo \
  --bam example/demo.bam \
  --ref example/genome.fa \
  --barcodes example/barcodes.txt \
  --vcf-dir example/vcfs \
  --out-dir example/output

Rscript scripts/combine_vartrix_matrices.R \
  --dir example/output \
  --vcf-dir example/vcf_sub \
  --barcodes example/barcodes.txt \
  --sample demo \
  --out-dir example/output \
  --plot-qc
```

---

## Notes

* Chromosomes 1–22 processed by default
* Missing files are skipped with warnings
* Output files are overwritten if present

---

## Use Cases

This pipeline supports:

* allele-specific CNV inference
* LOH detection
* tumor evolution analysis
* integration with CNVeil or similar frameworks

---

## Acknowledgments

* VarTrix

---

## Contact

Please open an issue for questions or bug reports.

# VarTrix Allele Matrix Pipeline

This script extracts **reference (REF)** and **alternative (ALT)** allele count matrices from a 10x Genomics BAM file using VarTrix.

It is designed for **single-cell DNA/RNA sequencing data** with **phased VCF files per chromosome**.

---

## Overview

Given:

* A 10x `possorted_bam.bam`
* A set of **phased VCF files** (per chromosome)
* A list of valid cell barcodes

This pipeline runs VarTrix across chromosomes and generates:

* REF allele count matrices
* ALT allele count matrices

These matrices can be used for downstream analysis such as:

* allele-specific copy number inference
* clonal structure analysis
* tools like Alleloscope / CNV methods

---

## Requirements

* VarTrix (tested with `vartrix_linux`)
* Samtools-indexed BAM file (`.bam + .bai`)
* Reference genome FASTA
* Phased VCF files (`.vcf.gz + .tbi`)
* Cell barcode list

---

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPO.git
cd YOUR_REPO
```

Make script executable:

```bash
chmod +x run_vartrix.sh
```

---

## Input Format

### BAM

10x Genomics BAM file:

```
possorted_bam.bam
```

---

### VCF files

Expected naming pattern:

```
<patient>_phased_chr1.vcf.gz
<patient>_phased_chr2.vcf.gz
...
```

Each VCF must be:

* bgzipped
* indexed (`.tbi`)

---

### Barcodes

A plain text file:

```
AAACCTGAGAAACCAT-1
AAACCTGAGAAACGAG-1
...
```

---

## Usage

```bash
bash run_vartrix.sh \
  --patient A \
  --bam breast_tissue_A_2k_possorted_bam.bam \
  --ref genome.fa \
  --barcodes uniq_cell.txt \
  --vcf-dir ./phased_vcfs \
  --out-dir ./vartrix_output \
  --vartrix /path/to/vartrix_linux
```

---

## Output

For each chromosome:

```
vartrix_output/
├── chr1_matrix_ref
├── chr1_matrix_alt
├── chr2_matrix_ref
├── chr2_matrix_alt
...
```

* `*_matrix_ref`: reference allele counts
* `*_matrix_alt`: alternative allele counts

---

## Notes

* Chromosomes 1–22 are processed by default
* Missing VCFs are skipped with a warning
* Existing output matrices are overwritten

---

## Tips

* Make sure BAM is indexed:

```bash
samtools index possorted_bam.bam
```

* Check VCF indexing:

```bash
tabix -p vcf file.vcf.gz
```

---

## Acknowledgments

This pipeline wraps functionality from:

* VarTrix (10x Genomics)

---

## Contact

For questions or issues, open a GitHub issue.

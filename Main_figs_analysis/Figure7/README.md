# Figure 7: Phylogenetic reconstruction and dosage sensitivity analysis

This directory contains the workflow used to reconstruct phylogenetic trees from total and haplotype-resolved copy number profiles and generate Figure 7.

The transcriptomic preprocessing workflow used for the dosage sensitivity analysis was adapted from the official wellDR-seq pipeline:

- [wellDR-seq GitHub repository](https://github.com/navinlabcode/wellDR-seq)

The original repository provides RNA preprocessing scripts, including the Perl FASTQ conversion script, STAR-based alignment workflow, and the chip barcode index file. 

## Prerequisites

Phylogenetic reconstruction was performed using MEDICC2:

- [MEDICC2]([https://github.com/cbg-ethz/MEDICC2](https://bitbucket.org/schwarzlab/medicc2/src/master/))

RNA preprocessing requires:

- Perl
- Java
- Trimmomatic
- STAR
- an hg19 FASTA file
- a corresponding GTF annotation
- the wellDR-seq barcode index file

Detailed installation instructions for MEDICC2 and the original wellDR-seq preprocessing pipeline are available from their respective repositories.

## Input

Required example files are provided under:

```text
example_data/figure7/
```

Main inputs include:

```text
example_data/figure7/
├── copy_number/
│   ├── input_medicc2.tsv
│   ├── subclone_totalcn.tsv
│   ├── total_cn_cluster_assignments.csv
│   ├── haplotype_cluster_assignments.tsv
│   └── total_cn_cluster_profiles.csv
├── transcriptomics/
│   ├── raw_fastq/
│   ├── com_index_wafer.txt
│   ├── cell_id_map.csv
│   ├── gene_count_matrix.tsv
│   └── gene_annotation.tsv
└── trees/
    ├── total_cn_support_tree.new
    └── haplotype_support_tree.new
```

### Barcode information

The wellDR-seq RNA reads require barcode-aware preprocessing before alignment.

The barcode index file:

```text
com_index_wafer.txt
```

contains the well or chip barcode definitions used by the wellDR-seq demultiplexing script. It is passed to:

```text
scDR-fq-convert-v4.pl
```

to separate RNA reads by cell.

The processed RNA cell barcodes are subsequently connected to the corresponding DNA cell identifiers using:

```text
cell_id_map.csv
```

This mapping is required to ensure that RNA expression profiles and DNA-derived subclone assignments refer to the same cells.

## Workflow

### 1. Preprocess the wellDR-seq RNA reads

```bash
bash scripts/04_preprocess_welldr_rna.sh
```

This step:

- demultiplexes the paired-end RNA FASTQ files using the wellDR-seq Perl script;
- uses `com_index_wafer.txt` to assign reads to individual cells;
- trims sequencing adapters;
- aligns the processed reads to hg19 using STAR;
- generates the gene-by-cell expression count matrix.

This preprocessing workflow was adapted from the `rna_fq_preprocess` workflow in the official wellDR-seq repository. The original repository includes:

```text
scDR-fq-convert-v4.pl
star-smart.sh
com_index_wafer.txt
```

Please refer to the original wellDR-seq repository for detailed information about the barcode layout and FASTQ conversion procedure.

### 2. Reconstruct phylogenetic trees

```bash
bash scripts/01_run_medicc2.sh
```

MEDICC2 is run separately using:

- haplotype-resolved copy number profiles;
- total copy number profiles.

This step generates the phylogenetic trees used in Figure 7a.

### 3. Compare total and haplotype-defined subclones

```bash
python scripts/02_phylogeny_and_cluster_correspondence.py
```

This step generates:

- Figure 7a: total and haplotype-resolved phylogenetic trees;
- Figure 7b: correspondence between total CN and haplotype-defined subclones.

### 4. Summarize cancer-gene copy number events

```bash
python scripts/03_cancer_gene_cn_events.py
```

This step generates:

- Figure 7c: cancer-gene total CN heatmap;
- Figure 7d: amplification, LOH, and CN-LOH frequencies across subclones.

### 5. Analyze dosage sensitivity

```bash
Rscript scripts/05_dosage_sensitivity.R
```

This step:

- matches RNA cell barcodes to DNA cell identifiers;
- aggregates RNA expression by inferred subclone;
- constructs pseudobulk expression profiles;
- maps CNVeil copy-number bins to genes;
- compares gene expression with integer copy number across subclones;
- generates the representative dosage-sensitive and dosage-insensitive gene panels shown in Figure 7e.

## Output

The workflow reproduces:

- **Figure 7a:** total and haplotype-specific phylogenetic trees;
- **Figure 7b:** total-CN to haplotype-cluster correspondence;
- **Figure 7c:** cancer-gene copy number heatmap;
- **Figure 7d:** cancer-gene amplification, LOH, and CN-LOH frequencies;
- **Figure 7e:** relationship between integer copy number and gene expression.

## Notes

The raw wellDR-seq sequencing data and the complete original preprocessing workflow are maintained by the wellDR-seq authors. This repository includes only the scripts and processed example inputs required to reproduce the CNVeil analyses and Figure 7.

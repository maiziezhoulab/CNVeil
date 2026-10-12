# Alleloscope: TN1 benchmark workflow

Three files cover Vartrix allele counting, combined-matrix preparation, and
Alleloscope execution. Benchmark evaluation and comparison plots live elsewhere.
The native Alleloscope diagnostic and lineage plots are retained.

## Requirements and inputs

- Bash, awk, sort, sed, and [Vartrix](https://github.com/10XGenomics/vartrix).
- R with [Alleloscope](https://github.com/seasoncloud/Alleloscope) and Matrix.
  The supplied cluster script used R 4.5.0; Alleloscope/Vartrix versions were not supplied.
- Indexed merged BAM with cell barcode CB tags and indexed matching genome FASTA.
- CHISEL `calls.tsv`: column 4 supplies cell IDs. The TN1 convention appends `-0`;
  verify these resulting barcodes exactly match the BAM CB tags.
- Autosomal phased VCFs named `TN1_phased_chr1.vcf` through
  `TN1_phased_chr22.vcf`, uncompressed. Variant calling/phasing happens upstream.
- Coverage CSV: first column is cell ID; other columns are numeric bin coverage.
  Bin names can be `chr1.0.500000` or `chr1-0-500000`.
  **Confirm that the original `cnv_matrix.csv` contains read coverage suitable
  for Alleloscope `raw_counts`, rather than inferred copy-number calls.**
- Normal-reference RDS: numeric bin-by-cell/bulk matrix for the same bins.
- GRCh37/hg19 chromosome-size table from the supplied Alleloscope installation
  (`sizes.cellranger-atac-hg19-1.2.0.txt`). BAM, FASTA, VCF, bins, and sizes must
  describe the same assembly. The preprocessing script adds `chr` to VCF contigs;
  it therefore assumes the BAM/FASTA use `chr` names.

## Run

Run from the directory containing these scripts. Output directories must be new;
their parent directories must exist for the Bash step.

```bash
VARTRIX=/path/to/vartrix bash prepare_vartrix.sh \
  TN1 /path/to/all_tumor_cells_norg.bam /path/to/genome.fa \
  /path/to/chisel/calls.tsv /path/to/phased_vcfs /path/to/vartrix_TN1

Rscript run.R TN1 /path/to/vartrix_TN1 \
  /path/to/cnv_matrix.csv /path/to/normalRD.rds \
  /path/to/sizes.cellranger-atac-hg19-1.2.0.txt /path/to/alleloscope_TN1
```

The BAM is supplied as an existing input; the commented RG-editing command in
the original script is not executed here. No BAM transformation is required by
these scripts. If your upstream workflow needs an edited BAM, prepare it there.

Optional settings:

| Setting | Default | Use |
|---|---|---|
| `ALLELOSCOPE_R_LIB` | Existing R library paths | Prepend a custom R package directory |
| `REFERENCE_REGION` | `8:5e+05` | Manual TN1 reference-region override |
| `COVERAGE_MULTIPLIER` | `2` | Scale genotype-value columns 1, 5, ..., 121 |
| `TRUST_NORMAL_BIN_ORDER` | Unset | Set to `1` only after verifying an unlabeled reference has identical bin order |
| `ALLELOSCOPE_CHECKPOINT` | Unset | Absolute path to a matching object after `Est_regions()` |

For cluster use, load your R module before execution, for example
`ml StdEnv/2023 r/4.5.0`. Original allocations were 10 CPUs/100 GB/3 days for
Alleloscope and 10 CPUs/100 GB/2 days for Vartrix. Pass site-specific account,
partition and QoS to `sbatch`; they are not embedded in the portable scripts.
For example, after configuring the module environment:

```bash
sbatch --cpus-per-task=10 --mem=100G --time=3-00:00:00 \
  --wrap='Rscript run.R TN1 /absolute/vartrix_TN1 /absolute/coverage.csv /absolute/normal.rds /absolute/sizes.tsv /absolute/new_results'
```

## Benchmark settings and changes

- Preserved: GRCh37, scDNAseq, `cell_filter=1000`, `SNP_filter=20`,
  VAF 0.1–0.9, segment `nSNP=2000`, `max_nSNP=30000`, `cont=TRUE`,
  `type="cellline"`, and lineage `nSNP=2000`, `nclust=10`.
- For TN1, the automatic reference selection yielded predominantly diploid
  calls, judged implausible by the benchmark author. The manual reference
  `8:5e+05` is retained; its diploid status has not been independently verified.
- The original ×2 adjustment to columns 1, 5, ..., 121 is retained exactly.
  Its biological rationale and column meanings still need confirmation.
  Set the multiplier to `1` to explicitly disable it. Another number/order of
  regions requires reviewing `coverage_columns`; the script does not guess.
- Corrected Vartrix labeling: `--out-matrix` writes ALT; `--ref-matrix` writes REF.
  The original commands reversed their filenames, and the combining function
  propagated that reversal. **The cleaned workflow can change inferred results.**
  Existing reversed matrices must not be used as inputs to this version.
- Preserved the effective original SNP inclusion rule: positive pooled REF
  coverage. The old function called this positive ALT coverage because the
  files were reversed. This is distinct from downstream VAF/SNP filtering.
- Cell columns now match by checked identifiers. SRR filename suffixes are
  stripped for matching, and the coverage-cell order is used consistently.
  Duplicate IDs after normalization or missing cells cause an error.
- Reference bins now match by identifiers. Unlabeled references need an explicit
  verified-order setting; labeled references with missing bins cause an error.
- Default execution now includes the previously commented object-building
  steps. Checkpoint mode is only appropriate for an object created with the same
  cells/bins and corrected allele counts. Do not resume a legacy reversed object
  to claim the allele correction was applied.
- Removed unused source/library calls, personal paths, debug prints, and
  commented alternatives. The original optional statistics histograms are omitted.

## Outputs and validation

Vartrix output contains the exact barcode list, transformed VCFs, and paired
chromosome matrices. Alleloscope output contains combined RDS inputs (before
coverage-cell subsetting), a pre-normal checkpoint for a fresh run, the object
before scaling, final `Obj_filtered.rds`, native plots, and `sessionInfo.txt`.

Bash syntax and a mocked 22-chromosome preprocessing run were checked during
organization. R and the biological inputs were unavailable in the editing
environment, so R syntax/package compatibility and end-to-end Alleloscope
execution remain unverified. Confirm the coverage input and scaling rationale
before treating this as a validated benchmark run.

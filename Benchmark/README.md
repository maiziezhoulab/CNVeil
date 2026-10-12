# SCOPE: hg19 workflow

```bash
Rscript run.R /path/to/hg19_individual_cell_bams /path/to/new_SCOPE_run
```

Requirements: your working SCOPE installation, WGSmapp,
BSgenome.Hsapiens.UCSC.hg19, and dependencies. Optional `SCOPE_R_LIB` prepends
a custom package library. All `.bam` files directly in the input directory
are selected; keep unrelated/merged BAMs outside it. Cell IDs are filename
prefixes before the first underscore and must be unique.

**Assembly change:** the original S0A script loaded hg38 and resumed from an
existing `qcObj.rds`. This version intentionally uses hg19, per the benchmark's
new shared requirement, and runs the previously commented preprocessing/QC
steps on actual hg19 BAMs. Do not use the original hg38 checkpoint or hg38
BAMs with these settings. This is a new analysis, not a coordinate conversion
or an exact reproduction of the original hg38 result.

| Preserved setting | Value |
|---|---|
| SCOPE resolution | 500 (500 kb) |
| Mapping quality threshold | 40 |
| Counting | paired-end |
| QC | `perform_qc()` installed defaults |
| Normal-reference selection | Gini ≤ 0.12 |
| Initial normalization | CODEX2 without latent factors |
| Ploidy | Per-cell `initialize_ploidy()` |
| SCOPE latent factors | K=1 |
| Iteration candidates | T=1:5 |
| Model selection | Maximum BIC |
| Cross-sample segmentation | `mode="integer"`, `max.ns=1` |

SCOPE's paired-end counting is retained even though FLCNA/HMMcopy used
single-end counting. That difference is explicit and should be reviewed if
counting settings are meant to be identical across this benchmark.
No new seed or grouped normalization is introduced.

Outputs: `qcObj.rds`, `segment_cs.rds`, `iCN.rds`, `input_files.txt`, and
`sessionInfo.txt`. `iCN` was already computed by the supplied code; it is now
also saved. The original job used R 4.2.1, 20 CPUs, 300 GB, and 7 days;
the segmentation loop remains sequential. R execution is unverified.

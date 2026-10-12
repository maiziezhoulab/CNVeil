# FLCNA: hg19 workflow

```bash
Rscript run.R /path/to/hg19_individual_cell_bams /path/to/new_FLCNA_TN1
```

Use the original working FLCNA installation with SCOPE, WGSmapp, and
BSgenome.Hsapiens.UCSC.hg19. Optional `FLCNA_R_LIB` prepends your R library.
The supplied FLCNA API uses `FLCNA(Y=...)` and `CNA.out(LRR=...)`; versions with
a different API stop with an explicit message rather than guessing matrix orientation.

The input directory must contain only the intended individual-cell hg19 BAMs.
All `.bam` files directly inside it are selected. Cell IDs are filename prefixes
before the first underscore; duplicates cause an error.
**The original filename pattern contained `hg38`, despite hg19 preprocessing.**
Confirm the BAMs were actually aligned to hg19. This script does not realign or
convert coordinates and does not infer the assembly from filenames.

| Preserved setting | Value |
|---|---|
| Reference | hg19 |
| SCOPE resolution | 500 (500 kb) |
| Mapping quality threshold | 40 |
| Counting | single-end |
| QC | mappability ≥ 0.9; GC range 20–80 |
| Candidate cluster counts | 3, 4, 5 |
| Lambda | 3 |
| CNA output | cutoff 0.35; L 100 |

Correction: `CNA.out()` now receives `QCobject$ref`, matching the bins retained
for normalization, instead of the unfiltered `ref_raw`. This can change output
when QC removes bins. Other supplied model choices and matrix orientation are
preserved. No additional sample QC or random seed is introduced.

Outputs: `output_FLCNA.rds`, `CNA.output.rds`, `retained_bins.rds`,
`input_files.txt`, and `sessionInfo.txt`.
Original cluster settings: R 4.2.1 via `ml GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1`,
50 CPUs, 100 GB, 7 days. The R code does not explicitly configure 50 workers.
R execution is unverified in the editing environment.

[FLCNA repository](https://github.com/FeifeiXiao-lab/FLCNA)

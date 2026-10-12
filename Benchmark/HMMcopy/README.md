# HMMcopy: hg19 workflow

```bash
Rscript run.R /path/to/hg19_individual_cell_bams \
  /path/to/original/correction.R /path/to/new_HMMcopy_TN1
```

Requirements: your working HMMcopy installation, SCOPE, WGSmapp,
BSgenome.Hsapiens.UCSC.hg19, data.table, and their dependencies.
Optional `HMMCOPY_R_LIB` prepends your package library.

**The original custom `correction.R` is required but was not provided.** It must
define `correctReadcount_github()`. The wrapper retains this function rather than
substituting the standard package correction, which might change results.

The input directory must contain only intended hg19 individual-cell BAMs;
all `.bam` files directly inside it are selected. IDs are filename prefixes
before the first underscore. No genome conversion is performed.

Preserved: hg19, SCOPE resolution 500 (500 kb), mapping quality 40,
single-end counting, GC divided by 100, exclusion of cells with zero total
coverage, custom read-count correction, and default `HMMsegment()` settings.

Cleaned: construct `chr/start/end/gc/map` by named fields, preserve matrix shape
even for one retained cell, retain cell names, and name the segment list.
`cell_names.rds` is now a character vector rather than a matrix/list built by
repeated `cbind`; adapt downstream readers if they expected that old shape.
The unused sample-QC calculation is omitted; it was not used to filter cells.

Outputs: `Y_raw.rds`, `Y.rds`, `ref_raw.rds`, `cell_names.rds`,
`segments_list_500k.rds`, `input_files.txt`, and `sessionInfo.txt`.
The duplicated `Y_raw_se` save is omitted. Original cluster allocation:
50 CPUs, 100 GB, 7 days; module `ml StdEnv/2023 r/4.5.0`.
The supplied segmentation loop remains sequential.
R/package/helper execution remains unverified without the original installation.

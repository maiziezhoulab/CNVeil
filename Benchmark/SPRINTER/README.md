# SPRINTER: hg19 workflow

```bash
bash run.sh /path/to/hg19_barcoded_tumor.bam /path/to/hg19_genome.fa \
  /path/to/new_SPRINTER_TN1
```

Requirements: Conda, a working CHISEL environment providing `chisel_rdr`, and
a working SPRINTER environment providing `sprinter`. The wrapper runs each via
`conda run` rather than requiring activation commands inside a batch shell.
Set `CHISEL_ENV` and `SPRINTER_ENV` to environment names or absolute environment
paths; defaults are `chisel` and `sprinter`.

Inputs are an indexed hg19 merged BAM with cell barcodes, an indexed matching
hg19 FASTA, and a new output directory. Both tools run inside that directory.
Actual coordinates and chromosome names must match; no genome conversion occurs.

| Preserved setting | Value |
|---|---|
| CHISEL RDR bin size | 50,000 bp |
| CHISEL minimum reads | 100,000 |
| CHISEL jobs | 8 |
| SPRINTER minimum reads | 100,000 |
| rtreads | 200 |
| cnreads | 1,000 |
| minnumcells | 15 |
| SPRINTER jobs | 8 |

All other settings retain installed defaults; no new seed is introduced.
The original cluster allocation was 20 CPUs, 100 GB, 2 days.

Outputs include `sprinter.input.tsv.gz`, both command logs, CHISEL RDR files,
and native SPRINTER outputs. Unlike the original script, this wrapper retains
intermediates instead of deleting `rdr/`, `baf/`, `combo/`, `calls/`, `clones/`,
and `plots/`. This requires additional disk space.
Bash syntax checked; execution and installed CLI compatibility remain unverified.

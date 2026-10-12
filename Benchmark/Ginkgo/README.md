# Ginkgo: hg19 workflow

```bash
bash run.sh /path/to/hg19_cell_beds /path/to/ginkgo /path/to/new_Ginkgo_TN1
```

Use your existing Ginkgo CLI installation with its compiled tools, R packages,
and `genomes/hg19/original/variable_500000_76_bowtie` reference bins.
The wrapper selects `*_hg19_rmdup.bed.gz`, matching the original selection.
BED coordinates must actually be hg19; filenames do not convert coordinates.

Preserved: `--genome hg19`, `--binning variable_500000_76_bowtie`, explicit cell
list, and defaults for other options. The commented alpha/bad-bin masking
options remain disabled. These are variable-width bins, not fixed 500 kb bins.

The native CLI writes outputs inside its input directory and can rewrite BED
files. This wrapper copies the selected BEDs into a fresh run directory first,
requiring extra disk space. Use paths without whitespace for the legacy CLI.
No new evaluation or plotting scripts are added.

Outputs: copied BEDs, `cell.list`, `ginkgo.log`, and native Ginkgo results,
including `SegCopy`. Original cluster allocation: 50 CPUs, 400 GB, 7 days.
The supplied CLI command did not explicitly configure 50 workers.
Bash syntax checked; native execution remains unverified.

[Ginkgo CLI source](https://github.com/robertaboukhalil/ginkgo/blob/master/cli/ginkgo.sh)

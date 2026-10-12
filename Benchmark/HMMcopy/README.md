# HMMcopy

Original tool: [HMMcopy](https://github.com/shahcompbio/HMMcopy)

## Run

Edit the paths at the top of the script, then run from this folder using your existing working installation.

```bash
Rscript run.R
```

For Slurm, use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` (add QoS if required). The original resource allocation and module command are retained in `submit.sh`.

## Inputs and settings

| Item | Value |
|---|---|
| Dependencies | HMMcopy, SCOPE, WGSmapp, hg19 BSgenome, data.table |
| Inputs | Individual-cell hg19 BAMs; original correction.R |
| Reference | hg19 |
| Coverage | 500 kb; MAPQ=40; single-end |
| Correction | correctReadcount_github(); GC/100 |
| Segmentation | HMMsegment() defaults; zero-count cells excluded |
| Outputs | Y_raw_se; Y_raw; Y; ref_raw; cell_names; segments_list_500k.rds |

Custom correction.R is required. Original reference column handling and intermediate file names are retained; zero-cell exclusion is simplified without changing its criterion.

These scripts receive a light cleanup of the supplied workflow. Actual data/package execution remains unverified.

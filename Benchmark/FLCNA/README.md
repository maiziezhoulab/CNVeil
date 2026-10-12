# FLCNA

Original tool: [FLCNA](https://github.com/FeifeiXiao-lab/FLCNA)

## Run

Edit the paths at the top of the script, then run from this folder using your existing working installation.

```bash
Rscript run.R
```

For Slurm, use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` (add QoS if required). The original resource allocation and module command are retained in `submit.sh`.

## Inputs and settings

| Item | Value |
|---|---|
| Dependencies | Original working FLCNA version; SCOPE, WGSmapp, hg19 BSgenome |
| Inputs | Individual-cell hg19 BAM directory |
| Reference | hg19 |
| Coverage | 500 kb; MAPQ=40; single-end |
| QC | mapp_thresh=0.9; gc_thresh=20–80 |
| Model | K=3,4,5; lambda=3; cutoff=0.35; L=100 |
| Outputs | output_FLCNA; CNA.output |

Original FLCNA API, orientation, and ref_raw output argument are retained. Check coordinate alignment if QC removes bins. The supplied hg38 filenames do not establish the actual BAM build; supply hg19 BAMs.

These scripts receive a light cleanup of the supplied workflow. Actual data/package execution remains unverified.

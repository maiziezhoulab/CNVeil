# AneuFinder

Original tool: [AneuFinder](https://github.com/ataudt/aneufinder)

## Run

Edit the paths at the top of the script, then run from this folder using your existing working installation.

```bash
Rscript run.R
```

For Slurm, use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` (add QoS if required). The original resource allocation and module command are retained in `submit.sh`.

## Inputs and settings

| Item | Value |
|---|---|
| Dependencies | AneuFinder; BSgenome.Hsapiens.UCSC.hg19 |
| Inputs | Merged tumor BAM for blacklist; individual-cell BAM directory |
| Reference | hg19; chr1–chr22 |
| Bins | 500 kb |
| Blacklist | ≤5th or ≥99.9th pooled count percentile |
| Model | GC correction; edivisive; numCPU=1; no breakpoint refinement |
| Outputs | Sample blacklist and native AneuFinder results |

Keep the pooled BAM separate from individual-cell inputs.

These scripts receive a light cleanup of the supplied workflow. Actual data/package execution remains unverified.

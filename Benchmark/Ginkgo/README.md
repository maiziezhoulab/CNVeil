# Ginkgo

Original tool: [Ginkgo](https://github.com/robertaboukhalil/ginkgo)

## Run

Edit the paths at the top of the script, then run from this folder using your existing working installation.

```bash
bash run.sh
```

For Slurm, use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` (add QoS if required). The original resource allocation and module command are retained in `submit.sh`.

## Inputs and settings

| Item | Value |
|---|---|
| Dependencies | Existing Ginkgo CLI and hg19 bin files |
| Inputs | Directory of *_hg19_rmdup.bed.gz files |
| Reference | hg19 |
| Binning | variable_500000_76_bowtie |
| Optional filters | Original defaults; maskbadbins not enabled |
| Outputs | cell.list and native results in the BED directory |

Original in-place execution is retained. Ginkgo can modify the BED directory; run on a working copy if needed.

These scripts receive a light cleanup of the supplied workflow. Actual data/package execution remains unverified.

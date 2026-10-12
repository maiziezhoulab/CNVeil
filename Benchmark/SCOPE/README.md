# SCOPE

Original tool: [SCOPE](https://github.com/rujinwang/SCOPE)

## Run

Edit the paths at the top of the script, then run from this folder using your existing working installation.

```bash
Rscript run.R
```

For Slurm, use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` (add QoS if required). The original resource allocation and module command are retained in `submit.sh`.

## Inputs and settings

| Item | Value |
|---|---|
| Dependencies | SCOPE, WGSmapp, hg19 BSgenome |
| Inputs | Existing hg19 qcObj.rds with Y, sampname, ref and QCmetric |
| Reference | hg19 checkpoint |
| Normal cells | Gini≤0.12 |
| Normalization | K=1; T=1:5; maximum BIC |
| Segmentation | integer; max.ns=1 |
| Outputs | segment_cs.rds |

Checkpoint-based execution is retained. The supplied S0A checkpoint was hg38; obtain a new hg19 checkpoint using hg19 BAMs. Changing the loaded BSgenome does not convert an existing checkpoint.

These scripts receive a light cleanup of the supplied workflow. Actual data/package execution remains unverified.

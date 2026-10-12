# SEACON

Original tool: [SEACON](https://github.com/NabaviLab/SEACON)

## Run

Edit the paths at the top of the script, then run from this folder using your existing working installation.

```bash
conda activate SEACON
bash run.sh
```

For Slurm, use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` (add QoS if required). The original resource allocation and module command are retained in `submit.sh`.

## Inputs and settings

| Item | Value |
|---|---|
| Dependencies | Existing SEACON environment and dependencies |
| Inputs | Cell BAM directory; hg19 FASTA/mappability; phased VCF; matched-normal BAM |
| Reference | hg19 |
| Read counts | 500 kb; P=12 |
| BAF blocks | 50 kb |
| Calling | tolerance=0.2; upper-filter=5; max-wgd=2 |
| Outputs | Native SEACON results in seacon_output |

Original prep_readcount → prep_baf → call sequence is retained. Confirm the matched-normal BAM is actually hg19 despite its original hg38 filename.

These scripts receive a light cleanup of the supplied workflow. Actual data/package execution remains unverified.

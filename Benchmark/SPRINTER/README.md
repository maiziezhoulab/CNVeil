# SPRINTER

Original tool: [SPRINTER](https://github.com/zaccaria-lab/SPRINTER)

## Run

Edit the paths at the top of the script, then run from this folder using your existing working installation.

```bash
bash run.sh
```

For Slurm, use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` (add QoS if required). The original resource allocation and module command are retained in `submit.sh`.

## Inputs and settings

| Item | Value |
|---|---|
| Dependencies | Conda environments providing CHISEL and SPRINTER |
| Inputs | Indexed hg19 barcoded BAM and matching FASTA |
| Reference | hg19 |
| RDR | 50 kb; minreads=100000; jobs=8 |
| SPRINTER | minreads=100000; rtreads=200; cnreads=1000 |
| Additional settings | minnumcells=15; jobs=8 |
| Outputs | sprinter.input.tsv.gz; native CHISEL RDR and SPRINTER results |

Original sequential Conda activation is retained. The broad directory deletion is omitted; intermediates remain in the selected output directory.

These scripts receive a light cleanup of the supplied workflow. Actual data/package execution remains unverified.

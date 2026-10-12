# CNRein

Original tool: [CNRein](https://github.com/elkebir-group/CNRein)

## Run

Edit the paths at the top of the script, then run from this folder using your existing working installation.

```bash
conda activate CNRein
bash run.sh
```

For Slurm, use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` (add QoS if required). The original resource allocation and module command are retained in `submit.sh`.

## Inputs and settings

| Item | Value |
|---|---|
| Dependencies | CNRein environment; samtools, bcftools, SHAPEIT4, GNU Parallel |
| Inputs | Raw BAM; original fix_bam_header.py; full hg19 reference directory; original CNRein repository |
| Reference | hg19; chr1–chr22 |
| Parallel jobs | 18 |
| SHAPEIT threads | 8 per job |
| Execution | Original cached stages, then full script.py |
| Outputs | Original intermediate directories and CNRein results |

fix_bam_header.py is still required. Original parallelism is preserved: SHAPEIT can use 144 threads despite the 30-CPU Slurm request. Do not reuse cached results from another input or assembly.

These scripts receive a light cleanup of the supplied workflow. Actual data/package execution remains unverified.

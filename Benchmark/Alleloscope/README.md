# Alleloscope

Original tool: [Alleloscope](https://github.com/seasoncloud/Alleloscope)

## Run

Edit the paths at the top of the script, then run from this folder using your existing working installation.

```bash
bash prepare_vartrix.sh
Rscript run.R
```

For Slurm, use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` (add QoS if required). The original resource allocation and module command are retained in `submit.sh`.

## Inputs and settings

| Item | Value |
|---|---|
| Dependencies | Alleloscope, DescTools, data.table, Matrix; Vartrix for preparation |
| Inputs | Combined allele RDS files, coverage CSV, normal RDS, hg19 sizes, and existing Obj_filtered checkpoint |
| Reference | GRCh37/hg19 |
| Manual region | 8:5e+05 |
| Genotype normalization | cellline; original columns 1,5,…,121 ×2 |
| Lineage | nSNP=2000; nclust=10 |
| Outputs | Obj_filtered_1; Obj_filtered.rds; native plots |

The supplied checkpoint-based execution is retained. The original Vartrix filenames are reversed relative to its ALT/REF output meaning; they are preserved here for fidelity, not endorsed as correct. Barcode order and normal-bin order must match the original inputs. Fresh object construction needs verified cell alignment.

These scripts receive a light cleanup of the supplied workflow. Actual data/package execution remains unverified.

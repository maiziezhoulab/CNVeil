# CHISEL

Original tool: [CHISEL](https://github.com/raphael-group/chisel)

## Run

Edit the paths at the top of the script, then run from this folder using your existing working installation.

```bash
bash phase.sh
bash run.sh
```

For Slurm, use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` (add QoS if required). The original resource allocation and module command are retained in `submit.sh`.

## Inputs and settings

| Item | Value |
|---|---|
| Dependencies | CHISEL, Eagle2, bcftools, bgzip, tabix |
| Inputs | Existing GATK VCF; hg19 phased genotype panels/map; rename.txt; indexed barcoded BAM and FASTA |
| Reference | hg19 |
| Eagle | 10 threads; uncompressed VCF |
| CHISEL | nonormal; seed=12; bin=500kb |
| Outputs | Phased chromosome VCFs; all_chromosomes.tsv; native CHISEL results |

Variant calling/BAM preparation remain upstream. The original TSV query is retained without new filtering; it assumes numeric chromosome names in Eagle output. Activate CHISEL before run.sh.

These scripts receive a light cleanup of the supplied workflow. Actual data/package execution remains unverified.

# SeCNV

Original tool: [SeCNV](https://github.com/deepomicslab/SeCNV)

## Run

Edit paths in `run.sh`. Activate your working SeCNV environment and use `bash run.sh`, or use `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` from this folder.

| Item | Value |
|---|---|
| Inputs | Directory of individual-cell BAMs; indexed hg19 FASTA |
| Software | Existing SeCNV Scripts directory and working environment |
| Reference assets | hg19 reference and matching mappability assets required by SeCNV |
| Bin size | 50,000 bp (50 kb), matching the supplied command |
| Other settings | Original installed-tool defaults |
| Output | Native SeCNV results in OUTPUT_DIR |

The original Scripts-directory execution and SeCNV.py command are retained. No additional model options or filtering are added. Execution remains unverified.

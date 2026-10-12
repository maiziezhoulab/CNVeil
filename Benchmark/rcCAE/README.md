# rcCAE

Original tool: [rcCAE](https://github.com/zhyu-lab/rccae)

## Run

Edit paths in `run.sh`, then use `bash run.sh` or `sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh` from this folder. Keep the original rcCAE installation and its `rccae` environment, compiled prepInput, and MATLAB Runtime v91.

| Item | Value |
|---|---|
| Inputs | Two BAM directories and their separate barcode lists |
| Reference | hg19 FASTA and hg19 mappability bigWig |
| Chromosomes | chr1–chr22 |
| Counting | 500 kb bins; mapping quality 0 |
| Count combination | First header; all columns from input 1; columns 4 onward from input 2 |
| Training | 50 epochs; batch size 64; lr 0.0001; latent dimension 3; seed 0 |
| SCHMM | MATLAB Runtime v91; final argument 10 |
| Outputs | readcounts1.txt, readcounts2.txt, readcounts.txt, lrc.txt and native results |

The supplied references were hg38; paths are changed to hg19 per the shared benchmark requirement. Both BAM inputs must actually use hg19. The original two-input combination is retained and assumes identical bin rows/order; the first header is not expanded. Logs now append both prepInput runs to log.txt. Execution remains unverified.

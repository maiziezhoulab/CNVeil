# Single-cell DNA CNV benchmark

This directory contains scripts for running 12 tools in a single-cell DNA copy-number benchmark. Each tool has its own folder with the run scripts, benchmark parameter settings, and Slurm submission script. Evaluation and comparison plotting are maintained separately.

## Supported tools

| Tool | Benchmark scripts | Original repository |
|---|---|---|
| Alleloscope | [Instructions](Alleloscope/README.md) | [GitHub](https://github.com/seasoncloud/Alleloscope) |
| AneuFinder | [Instructions](AneuFinder/README.md) | [GitHub](https://github.com/ataudt/aneufinder) |
| CHISEL | [Instructions](CHISEL/README.md) | [GitHub](https://github.com/raphael-group/chisel) |
| CNRein | [Instructions](CNRein/README.md) | [GitHub](https://github.com/elkebir-group/CNRein) |
| FLCNA | [Instructions](FLCNA/README.md) | [GitHub](https://github.com/FeifeiXiao-lab/FLCNA) |
| Ginkgo | [Instructions](Ginkgo/README.md) | [GitHub](https://github.com/robertaboukhalil/ginkgo) |
| HMMcopy | [Instructions](HMMcopy/README.md) | [GitHub](https://github.com/shahcompbio/HMMcopy) |
| rcCAE | [Instructions](rcCAE/README.md) | [GitHub](https://github.com/zhyu-lab/rccae) |
| SCOPE | [Instructions](SCOPE/README.md) | [GitHub](https://github.com/rujinwang/SCOPE) |
| SeCNV | [Instructions](SeCNV/README.md) | [GitHub](https://github.com/deepomicslab/SeCNV) |
| SEACON | [Instructions](SEACON/README.md) | [GitHub](https://github.com/NabaviLab/SEACON) |
| SPRINTER | [Instructions](SPRINTER/README.md) | [GitHub](https://github.com/zaccaria-lab/SPRINTER) |

## Inputs and reference genome

Input requirements differ by tool: some workflows use individual-cell BAMs, while others use merged barcoded BAMs, BED files, phased variants, or precomputed matrices. Check the tool-specific instructions before running. Alleloscope and SCOPE resume from precomputed objects by default.

## Running a tool

1. Install the tool and its dependencies using the linked repository instructions.
2. Open its folder and read the input requirements and benchmark settings.
3. Edit the data, reference, software, and output paths at the top of its scripts.
4. Activate the required environment and run the commands listed in its README.

For example, to run FLCNA:

```bash
cd FLCNA
Rscript run.R
```

For Slurm submission, use the provided wrapper from the tool's folder:

```bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION submit.sh
```

Adjust the module commands, resources, and QoS for your cluster. Tools with separate preprocessing scripts list their execution order in the corresponding README.

## Results

Each workflow saves the tool's native outputs to the location specified in its script. Output filenames and the parameters used in this benchmark are listed in the tool-specific README. Keep results from different datasets in separate output directories.

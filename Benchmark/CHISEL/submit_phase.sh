#!/bin/bash
#SBATCH --job-name=single_phase_1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=14-00:00:00
#SBATCH --output=12.out
#SBATCH --error=12.err
# Load your bcftools/htslib environment; supply account/partition on submission.
bash phase.sh

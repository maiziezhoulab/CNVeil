#!/bin/bash
#SBATCH --job-name=rccae_KTN
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --output=10.out
#SBATCH --error=10.err
# Supply your account and partition on submission.
bash run.sh

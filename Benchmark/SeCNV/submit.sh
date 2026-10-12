#!/bin/bash
#SBATCH --job-name=secnv_TN1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=500G
#SBATCH --time=7-00:00:00
#SBATCH --output=TN1_re.out
#SBATCH --error=TN1_re.err
# Supply your account and partition on submission.
ml GCC/10.2.0 SAMtools/1.12
ml Anaconda3/5.0.1
source activate SeCNV
bash run.sh

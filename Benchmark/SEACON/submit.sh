#!/bin/bash
#SBATCH --job-name=SEACON
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --output=10.out
#SBATCH --error=10.err
# Add your cluster account/partition/QoS options when submitting.
source ~/.bashrc
conda activate SEACON
bash run.sh
echo "SLURM_JOBID: $SLURM_JOBID"

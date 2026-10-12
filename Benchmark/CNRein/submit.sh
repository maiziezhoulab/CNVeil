#!/bin/bash
#SBATCH --job-name=CNRein
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=400G
#SBATCH --time=14-00:00:00
#SBATCH --output=10.out
#SBATCH --error=10.err
# Add your cluster account/partition/QoS options when submitting.
source ~/.bashrc
conda activate CNRein
bash run.sh
echo "SLURM_JOBID: $SLURM_JOBID"

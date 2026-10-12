#!/bin/bash
#SBATCH --job-name=SPRINTER
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --output=10.out
#SBATCH --error=10.err
# Add your cluster account/partition/QoS options when submitting.
# Conda environments are activated in run.sh.
bash run.sh
echo "SLURM_JOBID: $SLURM_JOBID"

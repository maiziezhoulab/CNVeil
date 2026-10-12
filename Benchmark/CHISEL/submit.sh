#!/bin/bash
#SBATCH --job-name=CHISEL
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --output=10.out
#SBATCH --error=10.err
# Add your cluster account/partition/QoS options when submitting.
# Activate CHISEL and load bcftools before submission.
bash run.sh
echo "SLURM_JOBID: $SLURM_JOBID"

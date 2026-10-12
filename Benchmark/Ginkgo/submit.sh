#!/bin/bash
#SBATCH --job-name=Ginkgo
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=400G
#SBATCH --time=7-00:00:00
#SBATCH --output=10.out
#SBATCH --error=10.err
# Add your cluster account/partition/QoS options when submitting.
# Load your existing Ginkgo environment.
bash run.sh
echo "SLURM_JOBID: $SLURM_JOBID"

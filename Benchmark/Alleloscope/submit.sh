#!/bin/bash
#SBATCH --job-name=Alleloscope
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=3-00:00:00
#SBATCH --output=10.out
#SBATCH --error=10.err
# Add your cluster account/partition/QoS options when submitting.
ml StdEnv/2023 r/4.5.0
Rscript --no-save run.R
echo "SLURM_JOBID: $SLURM_JOBID"

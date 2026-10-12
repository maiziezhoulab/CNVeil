#!/bin/bash
#SBATCH --job-name=FLCNA
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --output=10.out
#SBATCH --error=10.err
# Add your cluster account/partition/QoS options when submitting.
ml GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1
Rscript --no-save run.R
echo "SLURM_JOBID: $SLURM_JOBID"

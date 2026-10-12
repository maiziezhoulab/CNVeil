#!/bin/bash
#SBATCH --job-name=SCOPE
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=300G
#SBATCH --time=7-00:00:00
#SBATCH --output=10.out
#SBATCH --error=10.err
# Add your cluster account/partition/QoS options when submitting.
ml GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1
Rscript --no-save run.R
echo "SLURM_JOBID: $SLURM_JOBID"

#!/bin/bash
#SBATCH --job-name=AneuFinder
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --output=10.out
#SBATCH --error=10.err
# Add your cluster account/partition/QoS options when submitting.
module restore R_v421
Rscript --no-save run.R
echo "SLURM_JOBID: $SLURM_JOBID"

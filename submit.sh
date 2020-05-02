#!/bin/bash

#SBATCH --qos=short
#SBATCH --job-name=dev_run
#SBATCH --account=coen
#SBATCH --output=MCMCepi-%j-%N.out
#SBATCH --error=MCMCepi-%j-%N.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=END
#SBATCH --mail-user=opfeffer

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia/1.3.0
module load hpc/2015
julia -p 16 sample_inverse_mean_mi.jl

#!/usr/bin/bash --login

#SBATCH -n=1
#SBATCH -c=1
#SBATCH --threads-per-core=2
#SBATCH -m=1G
#SBATCH -t=00:01:00
#SBATCH --account=fl3
#SBATCH --partition=work

# Array test for SLURM

srun echo "Testing array, job $SLURM_ARRAY_TASK_ID"

#!/bin/bash
#SBATCH --job-name=torusefo
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=11:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/effobc_%A_%a.out
# Open-chain DMRG entropies, one line of effobc_tasks.txt per array task.
cd /scratch/lustre/jkudlerflam/torus_routeB/contrib/torus_routeB
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" effobc_tasks.txt)
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
python3 -u run_obc.py $LINE --chis 64 128 160

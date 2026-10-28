#!/bin/bash
#SBATCH --job-name=enceff
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-7
#SBATCH --output=logs/fes_%A_%a.out
# VUMPS finite-entanglement scaling, one line of eff_tasks.txt per array task.
cd /scratch/lustre/jkudlerflam/encoded_chain_numerics/encoded_chain_numerics
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" eff_tasks.txt)
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
python3 -u run_fes.py $LINE --maxit 600 --tol 1e-8

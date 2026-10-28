#!/bin/bash
#SBATCH --job-name=encfes
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=4-00:00:00
#SBATCH --array=0-20
#SBATCH --output=logs/fes_%A_%a.out
# VUMPS finite-entanglement scaling, one line of fes_tasks.txt per array task.
cd /scratch/lustre/jkudlerflam/encoded_chain_numerics/encoded_chain_numerics
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" fes_tasks.txt)
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
python3 -u run_fes.py $LINE --maxit 600 --tol 1e-8

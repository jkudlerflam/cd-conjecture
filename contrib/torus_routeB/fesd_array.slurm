#!/bin/bash
#SBATCH --job-name=torusfsd
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=11:00:00
#SBATCH --array=0-6
#SBATCH --output=logs/fesd_%A_%a.out
# VUMPS finite-entanglement scaling, one line of fesd_tasks.txt per array task.
cd /scratch/lustre/jkudlerflam/torus_routeB/contrib/torus_routeB
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" fesd_tasks.txt)
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
python3 -u run_vumps.py $LINE --chis 48 64 96 128 192 256 --maxit 400 --tol 1e-8

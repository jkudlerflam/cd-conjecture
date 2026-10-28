#!/bin/bash
#SBATCH --job-name=toruseff
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=11:00:00
#SBATCH --array=0-9
#SBATCH --output=logs/eff_%A_%a.out
# VUMPS finite-entanglement scaling, one line of eff_tasks.txt per array task.
cd /scratch/lustre/jkudlerflam/torus_routeB/contrib/torus_routeB
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" eff_tasks.txt)
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
python3 -u run_vumps.py $LINE --ncells 2 --chis 48 64 96 128 192 --maxit 400 --tol 1e-8

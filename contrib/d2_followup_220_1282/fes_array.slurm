#!/bin/bash
#SBATCH --job-name=d2fu_fes
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-15
#SBATCH --output=logs/fes_%A_%a.out
# Finite-entanglement scaling with VUMPS: one (model, cell, mode) per array task, see fes_tasks.txt.
# (Jobs 17511443 and 17511561.  Ising was run on a laptop instead of task 0.)
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
read MODEL N MODE <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" fes_tasks.txt)"
python3 run_fes.py "$MODEL" --n "$N" --mode "$MODE" --chis 8 12 16 24 32 48 64 96 128 \
    --maxit 3000 --tol 1e-8 --save_state --out "data/fes_${MODEL}_n${N}_${MODE}.json"

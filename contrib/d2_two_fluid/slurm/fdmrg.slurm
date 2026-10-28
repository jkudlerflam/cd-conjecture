#!/bin/bash
#SBATCH --job-name=tf_fd
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=03:00:00
#SBATCH --output=logs/fd_%A_%a.out
# Bulk finite-entanglement scaling on long open chains, one (name, spec, L) per task of file $1.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
read NAME SPEC L <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $1)"
python3 run_fdmrg.py "$SPEC" --L $L --chis 8 12 16 24 32 --win 40 80 --out data/fd_${NAME}_L${L}.json

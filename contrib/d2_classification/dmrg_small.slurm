#!/bin/bash
#SBATCH --job-name=d2dmrgS
#SBATCH --partition=typhon
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=1:30:00
#SBATCH --output=logs/dmrgS_%A_%a.out
# Small backfill-friendly DMRG tasks (chi <= 128), one line of tasks_dmrg_small.txt each.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
read NAME L <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" tasks_dmrg_small.txt)"
python3 dmrg.py "$NAME" --L "$L" --seed 3 --chi 32 64 128 --sweeps 2 --out "results/dmrg_${NAME}_L${L}_s3.json"

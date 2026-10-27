#!/bin/bash
#SBATCH --job-name=d2class
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=24
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --array=0-26
#SBATCH --output=logs/%A_%a.out
# Momentum-resolved ED of one chain per array task, L = 8..26, sectors in parallel.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" tasks_cluster.txt)
python3 run_point.py "$NAME" --procs 24 --sizes 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26

#!/bin/bash
#SBATCH --job-name=d2cls24
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=8:00:00
#SBATCH --output=logs/%A_%a.out
# Same as ed_array.slurm but L <= 24, sized to fit into partially allocated nodes.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" tasks_cluster.txt)
python3 run_point.py "$NAME" --procs 8 --sizes 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24

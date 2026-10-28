#!/bin/bash
#SBATCH --job-name=tf_tda
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --array=0-199
#SBATCH --output=logs/tda_%A_%a.out
# Two-body TDA screen of the exact-vacuum family: 200 tasks x 100 random points.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 screen_tda.py $((1000+SLURM_ARRAY_TASK_ID)) 100 data/tda_${SLURM_ARRAY_TASK_ID}.json

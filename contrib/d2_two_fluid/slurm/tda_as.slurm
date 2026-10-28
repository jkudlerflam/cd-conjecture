#!/bin/bash
#SBATCH --job-name=tf_tdas
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --array=0-39
#SBATCH --output=logs/tdas_%A_%a.out
# TDA screen of the antisymmetric-conversion slice gR = -gL e^{i delta}: 40 tasks x 60 points.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 screen_tda.py $((5000+SLURM_ARRAY_TASK_ID)) 60 data/tdas_${SLURM_ARRAY_TASK_ID}.json as

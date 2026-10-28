#!/bin/bash
#SBATCH --job-name=tf_sed
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/sed_%A_%a.out
# Exact ring check (screen_ed.py) of the TDA survivors at L = 12 and 14.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
LS=(12 14)
python3 screen_ed.py data/tda_random_survivors.json data/ed_survivors_L${LS[$SLURM_ARRAY_TASK_ID]}.json ${LS[$SLURM_ARRAY_TASK_ID]}

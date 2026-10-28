#!/bin/bash
#SBATCH --job-name=rA_sw3
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=0-01:30:00
#SBATCH --array=0-7
#SBATCH --output=logs/sw3_%A_%a.out
# Exact three-cell Schrieffer-Wolff for the C8 cell, one coupling per array task.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
LAMS=(0.0025 0.005 0.01 0.02 0.05 0.1 0.15 0.2)
L=${LAMS[$SLURM_ARRAY_TASK_ID]}
python3 sw_cluster.py --n3 --lams $L --out data/heff_c8_n3_lam${L}.npz

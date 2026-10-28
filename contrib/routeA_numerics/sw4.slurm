#!/bin/bash
#SBATCH --job-name=rA_sw4
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-12:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/sw4_%A_%a.out
# Exact four-cell full-chain low band versus the logical model (cluster-expansion check).
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
LAMS=(0.1 0.2)
python3 sw4_check.py ${LAMS[$SLURM_ARRAY_TASK_ID]}

#!/bin/bash
#SBATCH --job-name=tf_cal
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-5
#SBATCH --output=logs/cal_%A_%a.out
# iDMRG finite-entanglement calibration on models with known c.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
M=(xx xxz05 xxzm05 xxh xx2 xx2s)
python3 run_idmrg.py ${M[$SLURM_ARRAY_TASK_ID]} --chis 16 24 32 48 64 96 128 --steps 600 --tol 1e-10 \
    --out data/cal_${M[$SLURM_ARRAY_TASK_ID]}.jsonl

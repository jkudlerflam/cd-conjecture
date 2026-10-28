#!/bin/bash
#SBATCH --job-name=tf_cl
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --array=0-7
#SBATCH --output=logs/cl_%A_%a.out
# Product-state preemption check (classical.py) of the TDA survivors, 8 chunks.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 classical.py data/tda_random_survivors.json data/cl_random_survivors_${SLURM_ARRAY_TASK_ID}.json $SLURM_ARRAY_TASK_ID 8

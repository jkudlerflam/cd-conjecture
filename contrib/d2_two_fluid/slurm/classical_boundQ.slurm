#!/bin/bash
#SBATCH --job-name=tf_cl
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --array=0-19
#SBATCH --output=logs/cl_%A_%a.out
# Product-state preemption check (classical.py) of the TDA survivors, 20 chunks.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 classical.py data/tda_random_boundQ.json data/cl_random_boundQ_${SLURM_ARRAY_TASK_ID}.json $SLURM_ARRAY_TASK_ID 20

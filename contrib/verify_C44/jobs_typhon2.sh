#!/bin/bash
#SBATCH -J vC44g
#SBATCH -p typhon
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 03:30:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH --array=0-2
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd $SLURM_SUBMIT_DIR
models=(M0 M1 M2)
python3 -B run_c.py ${models[$SLURM_ARRAY_TASK_ID]} finite 192 48

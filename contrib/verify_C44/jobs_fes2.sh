#!/bin/bash
#SBATCH -J vC44w
#SBATCH -p typhon
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 02:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH --array=0-1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cd $SLURM_SUBMIT_DIR
models=(M1 M2)
python3 -B vumps_fes.py ${models[$SLURM_ARRAY_TASK_ID]} 96

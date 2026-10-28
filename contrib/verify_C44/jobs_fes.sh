#!/bin/bash
#SBATCH -J vC44v
#SBATCH -p typhon
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 03:30:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH --array=0-19
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cd $SLURM_SUBMIT_DIR
models=(M0 M1 M2 C0 C1)
chis=(64 80 96 128)
i=$SLURM_ARRAY_TASK_ID
python3 -B vumps_fes.py ${models[$((i / 4))]} ${chis[$((i % 4))]}

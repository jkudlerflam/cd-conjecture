#!/bin/bash
#SBATCH -J vC44f
#SBATCH -p typhon
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH -t 03:30:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH --array=0-9
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cd $SLURM_SUBMIT_DIR
models=(M0 M1 M2 C0 C1)
i=$SLURM_ARRAY_TASK_ID
m=${models[$((i / 2))]}
if [ $((i % 2)) -eq 0 ]; then python3 -B run_c.py $m finite 128 48; else python3 -B run_c.py $m finite 160 96; fi

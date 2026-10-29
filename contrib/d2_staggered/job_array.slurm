#!/bin/bash
#SBATCH --job-name=d2stag
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x_%A_%a.out
# one line of $TASKS (default tasks_main.txt) per array task; run from contrib/d2_staggered
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
TASKS=${TASKS:-tasks_main.txt}
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $TASKS)
python3 $LINE

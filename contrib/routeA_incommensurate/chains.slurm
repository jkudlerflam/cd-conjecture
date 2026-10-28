#!/bin/bash
#SBATCH --job-name=rAi_ch
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-20:00:00
#SBATCH --output=logs/ch_%A_%a.out
# Effective and calibration chains (chains.py), one line of $TASKS per array task.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
TASKS=${TASKS:-tasks_cal.txt}
LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $TASKS)"
echo "$LINE"
python3 -u chains.py $LINE

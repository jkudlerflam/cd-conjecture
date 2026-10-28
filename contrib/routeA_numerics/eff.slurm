#!/bin/bash
#SBATCH --job-name=rA_eff
#SBATCH --partition=apollo,typhon
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=0-23:00:00
#SBATCH --output=logs/eff_%A_%a.out
# Effective logical chains and calibrations; one line of $TASKS per array task.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
TASKS=${TASKS:-tasks_eff.txt}
LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $TASKS)"
echo "$LINE"
python3 effchain.py $LINE

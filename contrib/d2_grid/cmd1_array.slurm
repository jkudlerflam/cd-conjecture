#!/bin/bash
#SBATCH --job-name=d2grid_scr
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/cmd_%A_%a.out
# usage: sbatch --array=0-N cmd_array.slurm TASKFILE   (one shell command per line, run from contrib/d2_grid)
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $1)
echo "$LINE"
eval "$LINE"

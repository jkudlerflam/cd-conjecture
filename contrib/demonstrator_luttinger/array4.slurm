#!/bin/bash
#SBATCH --job-name=lutt4
#SBATCH --partition=typhon
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=23:50:00
#SBATCH --output=logs/%x_%a.out
# Typhon array over tasks_typhon4.txt (one line per task); run from the demonstrator_luttinger folder.
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" tasks_typhon4.txt)
python3 -u $LINE

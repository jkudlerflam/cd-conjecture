#!/bin/bash
#SBATCH -J fba_fam
#SBATCH -p apollo
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 24:00:00
#SBATCH --mem=120G
#SBATCH -o typhon/slurm_%j.out
#SBATCH -e typhon/slurm_%j.err
# Runs every line of the task file with 48 workers on one node.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
mkdir -p data/typhon
TASKS=${1:-typhon/tasks_stage1.txt}
cat "$TASKS" | xargs -P 32 -I{} bash -c '{} >> typhon/tasks_$SLURM_JOB_ID.log 2>&1'

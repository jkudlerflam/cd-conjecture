#!/bin/bash
#SBATCH --job-name=cdscan
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x_%A_%a.out
# usage: sbatch --array=0-N --mem=XXG --time=HH:MM:SS slurm_run.sh tasks.txt
# each line of tasks.txt is the argument list for run.py
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
cd /scratch/lustre/jkudlerflam/cdsearch
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$1")
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
/usr/bin/time -v python3 -u scan.py $LINE 2>&1 | grep -v "^\s*\(Page\|Socket\|Signals\|Swaps\|File\|Average\|Exit\|Minor\|Major\|Volun\|Invol\)"

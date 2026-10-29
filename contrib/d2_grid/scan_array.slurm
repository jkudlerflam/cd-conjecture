#!/bin/bash
#SBATCH --job-name=d2grid_scan
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/scan_%A_%a.out
# one line of scan_tasks.txt per array task; run from contrib/d2_grid
export PYTHONNOUSERSITE=1
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" scan_tasks.txt)
python3 scan.py $LINE --procs $SLURM_CPUS_PER_TASK

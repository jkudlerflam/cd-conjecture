#!/bin/bash
#SBATCH -J mcells_ed
#SBATCH -p typhon
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o logs/ed_%A_%a.log
#SBATCH --array=0-5
# Mobile-selector ED at L = 18 on Typhon (contrib/mobile_cells).  Tasks: points A, B x W = 10, 12, 14, S^z = 0.
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
cd /scratch/lustre/jkudlerflam/mobile_cells
mapfile -t T < tasks_ed.txt
python3 ed_one.py ${T[$SLURM_ARRAY_TASK_ID]}

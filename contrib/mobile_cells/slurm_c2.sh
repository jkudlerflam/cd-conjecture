#!/bin/bash
#SBATCH -J mcells_c2
#SBATCH -p typhon
#SBATCH -c 4
#SBATCH --mem=160G
#SBATCH -t 48:00:00
#SBATCH -o logs/c2_%A_%a.log
#SBATCH --array=0-2
# Two-flavour payload (c_p = 2) for contrib/mobile_cells: FES and ring ground states.
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
cd /scratch/lustre/jkudlerflam/mobile_cells
mapfile -t T < tasks_c2.txt
python3 ${T[$SLURM_ARRAY_TASK_ID]}

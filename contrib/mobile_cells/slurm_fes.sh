#!/bin/bash
#SBATCH -J mcells_fes
#SBATCH -p typhon
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 48:00:00
#SBATCH -o logs/fes_%A_%a.log
#SBATCH --array=0-5
# Finite-entanglement scaling runs for contrib/mobile_cells (VUMPS, numpy/scipy).
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
cd /scratch/lustre/jkudlerflam/mobile_cells
mapfile -t T < tasks_fes.txt
python3 run_fes.py ${T[$SLURM_ARRAY_TASK_ID]} --maxit 600

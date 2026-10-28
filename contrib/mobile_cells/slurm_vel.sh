#!/bin/bash
#SBATCH -J mcells_vel
#SBATCH -p typhon
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -t 24:00:00
#SBATCH -o logs/vel_%A_%a.log
#SBATCH --array=0-9
# Mode velocities and ground energies (z = 1, Casimir, uniqueness) at L = 16, 17, 18 for contrib/mobile_cells.
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
cd /scratch/lustre/jkudlerflam/mobile_cells
mapfile -t T < tasks_vel.txt
python3 ${T[$SLURM_ARRAY_TASK_ID]}

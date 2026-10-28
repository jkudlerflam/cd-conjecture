#!/bin/bash
#SBATCH --job-name=ulsSW
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=0-16:00:00
#SBATCH --array=0-19
#SBATCH --output=logs/uls_%A_%a.out
# Schur-Weyl spectra of the SU(d) ULS ring, one line of uls_tasks.txt per array task.
cd /scratch/lustre/jkudlerflam/refined_conjecture
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" uls_tasks.txt)
echo "task $SLURM_ARRAY_TASK_ID on $(hostname): $LINE"
python3 -u uls_schur_weyl.py $LINE

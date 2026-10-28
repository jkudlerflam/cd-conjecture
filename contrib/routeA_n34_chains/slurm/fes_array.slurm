#!/bin/bash
#SBATCH --job-name=rA34_fes
#SBATCH --partition=typhon
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=23:50:00
#SBATCH --output=logs/fes_%A_%a.out
# One line of slurm/fes_tasks.txt per array index (0-based).
cd /scratch/lustre/jkudlerflam/routeA_n34_chains
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" slurm/fes_tasks.txt)
echo "task $SLURM_ARRAY_TASK_ID: $LINE"
python3 -B run_fes34.py $LINE --nc 2 --mode damped

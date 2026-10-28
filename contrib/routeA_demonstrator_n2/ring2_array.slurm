#!/bin/bash
#SBATCH --job-name=rAn2_ring
#SBATCH --partition=typhon
#SBATCH --array=0-8
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/ring2_%a.out
cd /scratch/lustre/jkudlerflam/routeA_demo_n2
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" tasks_ring2.txt)
python3 ring2.py $LINE --out data/ring/ring2_${SLURM_ARRAY_TASK_ID}.json

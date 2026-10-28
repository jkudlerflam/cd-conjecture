#!/bin/bash
#SBATCH --job-name=rAdemo_ring
#SBATCH --partition=typhon
#SBATCH --array=0-17
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ring_%a.out
cd /scratch/lustre/jkudlerflam/routeA_demo
export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" tasks_ring.txt)
TAG=$(echo $LINE | tr -d ' -' )
python3 ring_ed.py $LINE --out data/ring/ring_${SLURM_ARRAY_TASK_ID}.json

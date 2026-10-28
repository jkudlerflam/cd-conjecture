#!/bin/bash
#SBATCH --job-name=rAdemo_fes
#SBATCH --partition=typhon
#SBATCH --array=0-12
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=logs/fes_%a.out
cd /scratch/lustre/jkudlerflam/routeA_demo
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" tasks_fes.txt)
NAME=$(echo $LINE | cut -d' ' -f1)
python3 run_fes.py $LINE --out data/fes/$NAME.json

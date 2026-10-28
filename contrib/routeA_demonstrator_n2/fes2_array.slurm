#!/bin/bash
#SBATCH --job-name=rAn2_fes
#SBATCH --partition=typhon
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/fes2_%a.out
cd /scratch/lustre/jkudlerflam/routeA_demo_n2
export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" tasks_fes2c.txt)
NAME=$(echo $LINE | cut -d' ' -f1)
python3 run_fes2.py $LINE --out data/fes/$NAME.json

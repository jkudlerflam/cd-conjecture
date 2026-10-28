#!/bin/bash
#SBATCH --job-name=rAn2_fesA
#SBATCH --partition=apollo
#SBATCH --array=0-2
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/fes2A_%a.out
cd /scratch/lustre/jkudlerflam/routeA_demo_n2
export OPENBLAS_NUM_THREADS=32 OMP_NUM_THREADS=32 MKL_NUM_THREADS=32
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" tasks_fes2c.txt)
NAME=$(echo $LINE | cut -d' ' -f1)
python3 run_fes2.py $LINE --out data/fes/${NAME}_A.json

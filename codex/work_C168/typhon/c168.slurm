#!/bin/bash
#SBATCH --job-name=c168
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c168/logs/c168_%j.out
cd /scratch/lustre/jkudlerflam/c168
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p logs out
timeout 2h python3 c168_job.py P 48 > logs/P.log 2>&1
timeout 1h python3 c168_job.py S1 48 > logs/S1.log 2>&1
timeout 4h python3 c168_job.py S2 48 > logs/S2.log 2>&1
python3 summarize.py out > summary.txt

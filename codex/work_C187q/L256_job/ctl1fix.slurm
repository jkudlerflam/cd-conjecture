#!/bin/bash
#SBATCH --job-name=c187c1
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c187b/logs/c187c1_%j.out
# Corrected c = 1 control (default t2) for the L = 256 stage; see L256_job/README.md amendment.
cd /scratch/lustre/jkudlerflam/c187b
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
timeout 350m python3 run.py cfg/ctl1fix_L256_N48.json out > logs/ctl1fix_L256_N48.log 2>&1

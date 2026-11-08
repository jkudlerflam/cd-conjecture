#!/bin/bash
#SBATCH --job-name=c152_scan
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c152/slurm_%j.out
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /scratch/lustre/jkudlerflam/c152
python3 typhon_scan.py 300 46

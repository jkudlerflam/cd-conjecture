#!/bin/bash
#SBATCH --job-name=C132-long
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=C132-long-%j.out
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
python3 run_checks.py --long

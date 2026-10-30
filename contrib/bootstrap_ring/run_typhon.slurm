#!/bin/bash
#SBATCH --job-name=bootring
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=60G
#SBATCH --time=23:00:00
#SBATCH --output=logs/%x-%j.out
# Usage: sbatch run_typhon.slurm tasks_file.txt
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=2
# Usage: sbatch --export=ALL,NPAR=24 run_typhon.slurm tasks.txt  (NPAR concurrent tasks)
mkdir -p out logs
python3 sweep.py "$1"

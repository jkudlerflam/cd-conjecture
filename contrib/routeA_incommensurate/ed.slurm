#!/bin/bash
#SBATCH --job-name=rAi_ed
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=0-12:00:00
#SBATCH --output=logs/ed_%j.out
# Two- and three-cell exact diagonalization of the full cell chain (ed_check.py).
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
python3 -u ed_check.py --target ${TARGET:-INC} --lams 0.0025 0.005 0.01 0.02 0.03 0.05 0.1 0.15 0.2 0.3 --lams3 0.005 0.01 0.02 0.05 0.1 0.2 0.3

#!/bin/bash
#SBATCH --job-name=c190a
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c190/logs/c190a_%j.out
# C190 (Claude Opus 5.5 agent, 2026-11-09), node A.  Full-qubit Z2 DMRG at the Fibonacci-cusp candidate bonds.
# 8 slots x 6 BLAS threads; field tuning first, then production (pipeline190.py, registered in RESULT.md).
cd /scratch/lustre/jkudlerflam/c18x/c190
mkdir -p logs outA
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 validate.py 10 > logs/validate_A.log 2>&1
python3 pipeline190.py --node A --slots 8 --threads 6 --out outA > logs/pipeline_A.log 2>&1
python3 analyze190.py outA > outA/summary.txt 2>&1

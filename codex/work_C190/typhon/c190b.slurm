#!/bin/bash
#SBATCH --job-name=c190b
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c190/logs/c190b_%j.out
# C190 (Claude Opus 5.5 agent, 2026-11-09), node B.  Full-qubit Z2 DMRG at the Fibonacci-cusp candidate bonds.
# 8 slots x 6 BLAS threads; field tuning first, then production (pipeline190.py, registered in RESULT.md).
cd /scratch/lustre/jkudlerflam/c18x/c190
mkdir -p logs outB
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 validate.py 10 > logs/validate_B.log 2>&1
python3 pipeline190.py --node B --slots 8 --threads 6 --out outB > logs/pipeline_B.log 2>&1
python3 analyze190.py outB > outB/summary.txt 2>&1

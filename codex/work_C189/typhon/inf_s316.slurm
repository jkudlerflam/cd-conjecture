#!/bin/bash
#SBATCH --job-name=c189_s316
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c189/logs/inf_s316_%j.out
# C189 infinite-chain run inf_s316 (Claude Opus 5.5 agent, 2026-11-09). idmrg.py via run_inf.py; see typhon/cfg/inf_s316.json.
cd /scratch/lustre/jkudlerflam/c18x/c189
mkdir -p logs out
export OMP_NUM_THREADS=10 OPENBLAS_NUM_THREADS=10 MKL_NUM_THREADS=10
python3 run_inf.py typhon/cfg/inf_s316.json out > logs/inf_s316.log 2>&1
python3 summarize_inf.py out > summary_inf_inf_s316.txt 2>&1

#!/bin/bash
#SBATCH --job-name=c189_s18
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c189/logs/inf_s18_%j.out
# C189 infinite-chain run inf_s18 (Claude Opus 5.5 agent, 2026-11-09). idmrg.py via run_inf.py; see typhon/cfg/inf_s18.json.
cd /scratch/lustre/jkudlerflam/c18x/c189
mkdir -p logs out
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
python3 run_inf.py typhon/cfg/inf_s18.json out > logs/inf_s18.log 2>&1
python3 summarize_inf.py out > summary_inf_inf_s18.txt 2>&1

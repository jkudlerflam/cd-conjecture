#!/bin/bash
#SBATCH --job-name=c189_ff2v
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c189/logs/inf_ff2v_%j.out
# C189 infinite-chain run inf_ff2v (Claude Opus 5.5 agent, 2026-11-09). idmrg.py via run_inf.py; see typhon/cfg/inf_ff2v.json.
cd /scratch/lustre/jkudlerflam/c18x/c189
mkdir -p logs out
export OMP_NUM_THREADS=10 OPENBLAS_NUM_THREADS=10 MKL_NUM_THREADS=10
python3 run_inf.py typhon/cfg/inf_ff2v.json out > logs/inf_ff2v.log 2>&1
python3 summarize_inf.py out > summary_inf_inf_ff2v.txt 2>&1

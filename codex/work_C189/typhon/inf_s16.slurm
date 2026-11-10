#!/bin/bash
#SBATCH --job-name=c189_s16
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c189/logs/inf_s16_%j.out
# C189 infinite-chain run inf_s16 (Claude Opus 5.5 agent, 2026-11-09). idmrg.py via run_inf.py; see typhon/cfg/inf_s16.json.
cd /scratch/lustre/jkudlerflam/c18x/c189
mkdir -p logs out
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
python3 run_inf.py typhon/cfg/inf_s16.json out > logs/inf_s16.log 2>&1
python3 summarize_inf.py out > summary_inf_inf_s16.txt 2>&1

#!/bin/bash
#SBATCH --job-name=c160
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c160/logs/c160_%j.out
# C160 packed job (Claude Opus agent, 2026-11-08).
#   D  36 open-chain DMRG runs (9 points of the K-stratum CSL family near h_F, L = 64, 128, 192, 256)  -> dmrg/*.json
#   X  12 classical texture tasks (4 random near-h_F bonds each, N = 2000, 16 starts)                   -> texture/*.jsonl
cd /scratch/lustre/jkudlerflam/c160
mkdir -p dmrg texture logs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
( for p in "h30_k040 -0.40 0.30" "h30_k070 -0.70 0.30" "h30_k085 -0.85 0.30" "h30_k095 -0.95 0.30" \
           "h30_k105 -1.05 0.30" "h30_k130 -1.30 0.30" "h00_k060 -0.60 0.00" "h45_k060 -0.60 0.45" "h45_k085 -0.85 0.45"; do
    for L in 256 192 128 64; do echo "D $p $L"; done
  done
  for s in $(seq 0 11); do echo "X $s"; done ) | \
  xargs -P 48 -L 1 sh -c 'if [ "$0" = D ]; then timeout 450m python3 dmrg_task.py $1 0.5 $2 $3 $4 dmrg > logs/D_$1_$4.log 2>&1; else timeout 450m python3 texture_task.py $1 4 2000 16 texture/tex_$1.jsonl > logs/X_$1.log 2>&1; fi'
python3 summarize.py > summary.txt 2>&1

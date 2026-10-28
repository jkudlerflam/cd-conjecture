#!/bin/bash
#SBATCH --job-name=rAi_ring
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=0-23:00:00
#SBATCH --output=logs/ring_%j.out
# Full-chain rings of L = 4 cells (D^4 = 16.8M), lowest levels per parity sector.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16
python3 -u bigchain.py ring --lam ${LAM:-0.1} --Ls 4 --nev 24 --ring_tol 1e-9 --out data/bigring4_INC_${LAM:-0.1}.json

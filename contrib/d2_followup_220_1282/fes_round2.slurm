#!/bin/bash
#SBATCH --job-name=d2fu_fes2
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/fes2_%j.out
# Round 2 of finite-entanglement scaling: continue the chi = 128 states of round 1 to chi = 192, 256.
# Four runs in parallel with 8 BLAS threads each.  Arguments: list of "MODEL N MODE" triples.
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
while [ $# -ge 3 ]; do
  MODEL=$1; N=$2; MODE=$3; shift 3
  cp "data/fes_${MODEL}_n${N}_${MODE}.json" "data/fes_${MODEL}_n${N}_${MODE}_r2.json"
  python3 run_fes.py "$MODEL" --n "$N" --mode "$MODE" --chis 128 192 256 --maxit 2000 --tol 1e-8 --save_state \
      --load "data/fes_${MODEL}_n${N}_${MODE}_chi128.npz" --out "data/fes_${MODEL}_n${N}_${MODE}_r2.json" \
      > "logs/fes2_${MODEL}_n${N}.log" 2>&1 &
done
wait

#!/bin/bash
# C190g gate G4 benchmark, run as an overlap step inside the tests allocation (Claude Opus 5.5 agent, 2026-11-09).
# Usage: bench_step.sh <GPU UUID>
cd /scratch/lustre/jkudlerflam/c18x/c190g
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export CUDA_VISIBLE_DEVICES=$1 OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
nvidia-smi -L
python3 bench.py cupy --L 256 --chis 256,384,512,768 --nb 3 --nmv 80 --cpu 1 --out logs/bench_L256.json > logs/bench_L256.log 2>&1
python3 bench.py cupy --L 256 --chis 768 --nb 2 --nmv 80 --svd-dev --out logs/bench_L256_svddev.json > logs/bench_L256_svddev.log 2>&1
python3 bench.py cupy --L 256 --chis 32,64,128,192 --nb 3 --nmv 80 --cpu 1 --out logs/bench_L256_small.json > logs/bench_L256_small.log 2>&1
echo bench_step done

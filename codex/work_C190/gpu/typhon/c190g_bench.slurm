#!/bin/bash
#SBATCH --job-name=c190g_bench
#SBATCH --partition=apollo
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=05:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c190g/logs/c190g_bench_%j.out
# C190g gate G4 (Claude Opus 5.5 agent, 2026-11-09): L = 256 complex (p153c) timing per matvec and per bond on the
# GPU at chi = 256 to 768, against the original CPU code (6 BLAS threads, as in production) on the same inputs.
cd /scratch/lustre/jkudlerflam/c18x/c190g
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
nvidia-smi
(while true; do nvidia-smi --query-gpu=name,memory.used,memory.total,utilization.gpu --format=csv,noheader >> logs/bench_smi.log; sleep 30; done) &
python3 bench.py cupy --L 256 --chis 256,384,512,768 --nb 3 --nmv 80 --cpu 1 --out logs/bench_L256.json > logs/bench_L256.log 2>&1
python3 bench.py cupy --L 256 --chis 768 --nb 2 --nmv 80 --svd-dev --out logs/bench_L256_svddev.json > logs/bench_L256_svddev.log 2>&1
python3 bench.py cupy --L 256 --chis 32,64,128,192 --nb 3 --nmv 80 --cpu 1 --out logs/bench_L256_small.json > logs/bench_L256_small.log 2>&1
kill %1

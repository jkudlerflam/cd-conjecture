#!/bin/bash
#SBATCH --job-name=c190g_g2
#SBATCH --partition=apollo
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c190g/logs/c190g_g2_%j.out
# C190g gate G2 (Claude Opus 5.5 agent, 2026-11-09): cupy backend against ED at L = 10, SVD on host and on device.
cd /scratch/lustre/jkudlerflam/c18x/c190g
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
nvidia-smi
python3 -c "import cupy; print(cupy.__version__, cupy.cuda.runtime.runtimeGetVersion(), cupy.cuda.runtime.driverGetVersion())"
python3 validate_xp.py cupy 10 --sweep32 > logs/G2_cupy_svdhost.log 2>&1
python3 validate_xp.py cupy 10 --svd-dev > logs/G2_cupy_svddev.log 2>&1

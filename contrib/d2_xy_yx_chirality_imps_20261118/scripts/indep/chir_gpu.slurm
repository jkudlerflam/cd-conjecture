#!/bin/bash
#SBATCH --partition=apollo
#SBATCH --qos=gpu-medium
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c190g_chir/logs/%x_%j.out
# Independent-code check of the four-term chain (2026-11-24): one run190g run on one GPU, resumable.
# usage: sbatch -J <tag> chir_gpu.slurm <tag>
cd /scratch/lustre/jkudlerflam/c18x/c190g_chir
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export C190G_BACKEND=cupy C190G_GPU_MIN_CHI=128 C190G_SVD_DEV=1 C190G_PROFILE=0 OMP_NUM_THREADS=6
echo "GPU $CUDA_VISIBLE_DEVICES on $(hostname) $(date)"
if ! timeout 300 python3 cupy_smoke.py; then
  export CUPY_CACHE_DIR=/scratch/lustre/jkudlerflam/c18x/c190g_chir/cupy_cache_$(hostname -s); mkdir -p "$CUPY_CACHE_DIR"
  timeout 300 python3 cupy_smoke.py || { echo "SMOKE FAILED"; exit 3; }
fi
python3 run_chir.py cfg/$1.json outC
echo "rc=$? $(date)"

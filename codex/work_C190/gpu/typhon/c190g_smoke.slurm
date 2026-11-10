#!/bin/bash
#SBATCH --job-name=c190g_smoke
#SBATCH --partition=apollo
#SBATCH --qos=gpu-medium
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c190g/logs/c190g_smoke_%j.out
cd /scratch/lustre/jkudlerflam/c18x/c190g
echo start $(date) $(hostname) $CUDA_VISIBLE_DEVICES
timeout 60 nvidia-smi -L; echo smi rc=$?
timeout 120 python3 -c "import numpy; print('numpy ok')"; echo rc=$?
timeout 120 python3 -c "import cupy; print('cupy import ok')"; echo rc=$?
timeout 120 python3 -c "import cupy; print(cupy.cuda.runtime.getDeviceCount())"; echo rc=$?
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
timeout 120 python3 -c "import cupy.cuda.compiler as c; print('nvrtc', c._get_nvrtc_version())"; echo rc=$?
timeout 300 python3 cupy_smoke.py; echo rc=$?
echo end $(date)

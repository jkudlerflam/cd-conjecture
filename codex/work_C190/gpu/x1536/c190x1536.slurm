#!/bin/bash
#SBATCH --job-name=c190x1536
#SBATCH --partition=apollo
#SBATCH --qos=gpu-short
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c19x/c190x1536/logs/c190x1536_%j.out
# chi = 1536 stage at L = 256 (C190g R13, coordinator, 2026-11-10). Two runs sharing one GPU (offloaded storage keeps
# each well under 20 GB of device memory), C195's offloaded solver.
# A resubmission resumes each run from out/ckpt and skips finished ones.
cd /scratch/lustre/jkudlerflam/c19x/c190x1536
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export C190G_BACKEND=cupy C190G_GPU_MIN_CHI=128 C190G_SVD_DEV=1 C190G_PROFILE=0 C195_OFFLOAD=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
echo "GPUs $CUDA_VISIBLE_DEVICES on $(hostname) $(date)"
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader
( python3 run195.py cfg/p153c_L256_P1_rand1.json out > logs/p153c_L256_P1_rand1.log 2>&1 ) &
( python3 run195.py cfg/sym40_n316_L256_P1_rand1.json out > logs/sym40_n316_L256_P1_rand1.log 2>&1 ) &
wait
echo "done $(date)"

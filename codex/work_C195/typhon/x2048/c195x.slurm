#!/bin/bash
#SBATCH --partition=apollo
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=260G
#SBATCH --output=/scratch/lustre/jkudlerflam/c19x/c195x2048/logs/%x_%j.out
# C195x (coordinator, 2026-11-11): chi = 2048 stage at L = 512, one run, one GPU, C195's offloaded solver.
# Usage: sbatch --job-name=c195x_<key> --qos=<qos> --time=<t> c195x.slurm <tag>. A resubmission resumes from out2048/ckpt.
cd /scratch/lustre/jkudlerflam/c19x/c195x2048
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export C190G_BACKEND=cupy C190G_GPU_MIN_CHI=128 C190G_SVD_DEV=1 C190G_PROFILE=0 C195_OFFLOAD=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
echo "GPU $CUDA_VISIBLE_DEVICES on $(hostname) $(date) job $SLURM_JOB_ID"
(while true; do echo "$(date +%FT%T) $(nvidia-smi --query-gpu=memory.used --format=csv,noheader) $(free -g | sed -n 2p)" >> logs/mem_$SLURM_JOB_ID.log; sleep 300; done) &
M=$!
python3 run195.py cfg/$1.json out2048 > logs/$1.log 2>&1
echo "rc=$? $(date)"
kill $M

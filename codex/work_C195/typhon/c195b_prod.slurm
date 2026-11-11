#!/bin/bash
#SBATCH --job-name=c195b
#SBATCH --partition=apollo
#SBATCH --qos=gpu-short
#SBATCH --gres=gpu:2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=320G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c19x/c195b/logs/c195b_%j.out
# C195b (coordinator, 2026-11-10): sym40_n316 at L = 512, both sectors, one run per GPU, C195's offloaded solver.
# A resubmission resumes every unfinished run from out512/ckpt and skips finished ones.
cd /scratch/lustre/jkudlerflam/c19x/c195b
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export C190G_BACKEND=cupy C190G_GPU_MIN_CHI=128 C190G_SVD_DEV=1 C190G_PROFILE=0 C195_OFFLOAD=1
echo "GPUs $CUDA_VISIBLE_DEVICES on $(hostname) $(date) job $SLURM_JOB_ID"
python3 pipeline195b.py --out out512 --threads 4 > logs/pipeline_$SLURM_JOB_ID.log 2>&1
echo "pipeline rc=$? $(date)"

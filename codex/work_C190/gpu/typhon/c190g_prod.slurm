#!/bin/bash
#SBATCH --job-name=c190g_prod
#SBATCH --partition=apollo
#SBATCH --qos=gpu-medium
#SBATCH --gres=gpu:6
#SBATCH --nodes=1
#SBATCH --nodelist=apollo02
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=320G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c190g/logs/c190g_prod_%j.out
# C190g production (Claude Opus 5.5 agent, 2026-11-09).  One run per GPU (pipeline190g.py), 4 threads per run.
# Output outG/ (same layout as the CPU outA/: <tag>.json, logs/<tag>.log, cfg/, ckpt/).  A resubmission resumes
# every unfinished run from its last checkpoint and skips finished ones.
cd /scratch/lustre/jkudlerflam/c18x/c190g
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export C190G_BACKEND=cupy C190G_GPU_MIN_CHI=128 C190G_SVD_DEV=1 C190G_PROFILE=0
echo "GPUs $CUDA_VISIBLE_DEVICES on $(hostname) $(date)"
nvidia-smi -L
(while true; do nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader >> logs/prod_smi_$SLURM_JOB_ID.log; sleep 300; done) &
SMI=$!
python3 pipeline190g.py --out outG --threads 4 > logs/pipeline_G_$SLURM_JOB_ID.log 2>&1
kill $SMI

#!/bin/bash
#SBATCH --job-name=c195_prod
#SBATCH --partition=apollo
#SBATCH --gres=gpu:2
#SBATCH --nodes=1
#SBATCH --nodelist=apollo02
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=320G
#SBATCH --output=/scratch/lustre/jkudlerflam/c19x/c195/logs/c195_prod_%j.out
# C195 production (Claude Opus 5.5 agent, 2026-11-10).  p153c at L = 512, P = 0 and 1, one run per GPU
# (pipeline195.py, run195.py with the offloaded storage), 4 threads per run, up to 160 GB of host memory per run.
# QOS and time limit are given on the sbatch command line (RESULT.md).  A resubmission of this script resumes every
# unfinished run from its checkpoint (out512/ckpt) and skips finished ones.
cd /scratch/lustre/jkudlerflam/c19x/c195
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export C190G_BACKEND=cupy C190G_GPU_MIN_CHI=128 C190G_SVD_DEV=1 C190G_PROFILE=0 C195_OFFLOAD=1
echo "GPUs $CUDA_VISIBLE_DEVICES on $(hostname) $(date) job $SLURM_JOB_ID"
nvidia-smi --query-gpu=index,name,memory.used,memory.total --format=csv,noheader
(while true; do echo "$(date +%FT%T)" >> logs/prod_smi_$SLURM_JOB_ID.log; nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader >> logs/prod_smi_$SLURM_JOB_ID.log; free -g | sed -n 2p >> logs/prod_smi_$SLURM_JOB_ID.log; sleep 300; done) &
SMI=$!
python3 pipeline195.py --out out512 --threads 4 > logs/pipeline_$SLURM_JOB_ID.log 2>&1
echo "pipeline rc=$? $(date)"
kill $SMI

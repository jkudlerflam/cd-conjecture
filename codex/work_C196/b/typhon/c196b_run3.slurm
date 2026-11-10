#!/bin/bash
#SBATCH --job-name=c196b
#SBATCH --partition=apollo
#SBATCH --qos=gpu-medium
#SBATCH --gres=gpu:3
#SBATCH --nodes=1
#SBATCH --nodelist=apollo02
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=190G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c19x/c196b/logs/c196b_%j.out
# C196b (Claude Opus 5.5 agent, 2026-11-10). Up to three runs of run196b.py, one per GPU, 4 threads each.
# Usage: sbatch --job-name=c196b_<x> c196b_run3.slurm <outdir> <cfg1> [<cfg2> [<cfg3>]]  (cfg names in cfg/).
# Resubmitting the same command resumes every run from its checkpoint and skips finished chi.
cd /scratch/lustre/jkudlerflam/c19x/c196b
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export C196B_BACKEND=cupy C196B_GPU_MIN_CHI=128 C196B_SVD_DEV=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
OUT=$1; shift
mkdir -p $OUT
echo "GPUs $CUDA_VISIBLE_DEVICES on $(hostname) $(date) cfgs $*"
IFS=',' read -ra G <<< "$CUDA_VISIBLE_DEVICES"
timeout 300 python3 cupy_smoke.py || { echo "SMOKE FAILED"; exit 3; }
(while true; do nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader >> logs/smi_$SLURM_JOB_ID.log; sleep 600; done) &
SMI=$!
i=0
for c in "$@"; do
  CUDA_VISIBLE_DEVICES=${G[$i]} python3 -E run196b.py cfg/$c.json $OUT > logs/${c}_$SLURM_JOB_ID.log 2>&1 &
  i=$((i+1))
done
wait $(jobs -p | grep -v "^$SMI$")
echo "done $(date)"
kill $SMI

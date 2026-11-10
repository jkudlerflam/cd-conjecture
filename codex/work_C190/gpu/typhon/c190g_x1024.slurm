#!/bin/bash
#SBATCH --job-name=c190g_x1024
#SBATCH --partition=h200
#SBATCH --qos=gpu-h200
#SBATCH --gres=gpu:3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c190g/logs/c190g_x1024_%j.out
# C190g, chi = 1024 stage at L = 256 (coordinator, Claude Fable 5, 2026-11-10).  Three runs at a time, one per GPU
# (pipeline_x1024.py), 4 threads per run.  Output outX1024/ (same layout as outG/).  Nothing in outG/ is written.
# A resubmission resumes every unfinished run from its checkpoint in outX1024/ckpt and skips finished ones.
cd /scratch/lustre/jkudlerflam/c18x/c190g
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export C190G_BACKEND=cupy C190G_GPU_MIN_CHI=128 C190G_SVD_DEV=1 C190G_PROFILE=0
echo "GPUs $CUDA_VISIBLE_DEVICES on $(hostname) $(date)"
which python3; python3 -V
timeout 60 nvidia-smi --query-gpu=index,name,memory.total,memory.used --format=csv,noheader
# The cupy path is untested on the h200 nodes, so test it before production (first with the default kernel
# cache, then with a cache folder of this node).
ok=0
if timeout 300 python3 cupy_smoke.py; then ok=1; else
  echo "smoke test failed with the default kernel cache (rc=$?), retrying with a cache folder of this node"
  export CUPY_CACHE_DIR=/scratch/lustre/jkudlerflam/c18x/c190g/cupy_cache_$(hostname -s)
  mkdir -p "$CUPY_CACHE_DIR"
  if timeout 300 python3 cupy_smoke.py; then ok=1; fi
fi
if [ "$ok" != 1 ]; then echo "SMOKE FAILED on $(hostname), no production started"; exit 3; fi
echo "smoke ok $(date)"
(while true; do nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader >> logs/x1024_smi_$SLURM_JOB_ID.log; sleep 300; done) &
SMI=$!
python3 pipeline_x1024.py --src outG --out outX1024 --threads 4 > logs/pipeline_X_$SLURM_JOB_ID.log 2>&1
echo "pipeline rc=$? $(date)"
kill $SMI

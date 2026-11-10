#!/bin/bash
#SBATCH --job-name=c190g_tests
#SBATCH --partition=apollo
#SBATCH --qos=gpu-medium
#SBATCH --gres=gpu:4
#SBATCH --nodes=1
#SBATCH --nodelist=apollo02
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c190g/logs/c190g_tests_%j.out
# C190g gates G2, G3, G4 in one allocation (Claude Opus 5.5 agent, 2026-11-09).  The gpu-short QOS allows two jobs
# per user and both are in use, so the tests share one 4-GPU job, one GPU each.
cd /scratch/lustre/jkudlerflam/c18x/c190g
# cupy-cuda12x needs CUDA headers for NVRTC; the only CUDA 12 headers on the system are the nvidia-cuda-runtime-cu12
# 12.6 wheel files in ~/.local (python 3.13), linked as cudahdr/include.  apollo01 hangs on the cupy import, so use apollo02.
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
IFS=',' read -ra G <<< "$CUDA_VISIBLE_DEVICES"
echo "GPUs: $CUDA_VISIBLE_DEVICES on $(hostname)"
nvidia-smi
python3 -c "import cupy; print(cupy.__version__, cupy.cuda.runtime.runtimeGetVersion(), cupy.cuda.runtime.driverGetVersion())"
(while true; do nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu --format=csv,noheader >> logs/tests_smi.log; sleep 60; done) &
SMI=$!
( export CUDA_VISIBLE_DEVICES=${G[0]} OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
  python3 validate_xp.py cupy 10 --sweep32 > logs/G2_cupy_svdhost.log 2>&1
  python3 validate_xp.py cupy 10 --svd-dev > logs/G2_cupy_svddev.log 2>&1 ) &
P1=$!
( export CUDA_VISIBLE_DEVICES=${G[1]} OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
  python3 bench.py cupy --L 256 --chis 256,384,512,768 --nb 3 --nmv 80 --cpu 1 --out logs/bench_L256.json > logs/bench_L256.log 2>&1
  python3 bench.py cupy --L 256 --chis 768 --nb 2 --nmv 80 --svd-dev --out logs/bench_L256_svddev.json > logs/bench_L256_svddev.log 2>&1
  python3 bench.py cupy --L 256 --chis 32,64,128,192 --nb 3 --nmv 80 --cpu 1 --out logs/bench_L256_small.json > logs/bench_L256_small.log 2>&1 ) &
P2=$!
mkdir -p g3out
( export CUDA_VISIBLE_DEVICES=${G[2]} OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 C190G_BACKEND=cupy C190G_PROFILE=1
  python3 run190g.py cfg_g3/p153c_L64_P0.json g3out > logs/G3_P0.log 2>&1 ) &
P3=$!
( export CUDA_VISIBLE_DEVICES=${G[3]} OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 C190G_BACKEND=cupy C190G_PROFILE=0
  python3 run190g.py cfg_g3/p153c_L64_P1.json g3out > logs/G3_P1.log 2>&1 ) &
P4=$!
wait $P1 $P2 $P3 $P4
kill $SMI

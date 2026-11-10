#!/bin/bash
#SBATCH --job-name=c195_gate
#SBATCH --partition=apollo
#SBATCH --qos=gpu-medium
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --nodelist=apollo02
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=05:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c19x/c195/logs/c195_gate_%j.out
# C195 gate and benchmark on one A100 (Claude Opus 5.5 agent, 2026-11-10).
# (1) smoke test of bench195.py at L = 64, (2) gate: p153c L = 128 short schedule with the unmodified C190g code
# (ref/run190g.py) and with the offloaded C195 code (run195.py), plus a C195 run stopped after 8 sweeps and resumed,
# three processes on the one GPU, (3) benchmark at L = 512, chi = 768, 1024, 1536 (alone on the GPU).
cd /scratch/lustre/jkudlerflam/c19x/c195
export CUDA_PATH=/scratch/lustre/jkudlerflam/c18x/c190g/cudahdr
export C190G_BACKEND=cupy C190G_GPU_MIN_CHI=128 C190G_SVD_DEV=1 C190G_PROFILE=0
echo "GPU $CUDA_VISIBLE_DEVICES on $(hostname) $(date)"
nvidia-smi --query-gpu=index,name,memory.used,memory.total --format=csv,noheader
(while true; do nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader >> logs/gate_smi_$SLURM_JOB_ID.log; free -g | sed -n 2p >> logs/gate_smi_$SLURM_JOB_ID.log; sleep 60; done) &
SMI=$!
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
timeout 900 python3 bench195.py --L 64 --chis 64 --nb 1 --nmv 6 > logs/smoke_bench_$SLURM_JOB_ID.log 2>&1
echo "smoke rc=$? $(date)"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p gate/ref gate/off gate/res
( python3 ref/run190g.py cfg/gate_p153c_L128_P0.json gate/ref > gate/ref.log 2>&1 ) &
P1=$!
( C195_OFFLOAD=1 python3 run195.py cfg/gate_p153c_L128_P0.json gate/off > gate/off.log 2>&1 ) &
P2=$!
( C195_OFFLOAD=1 C195_EXIT_AFTER=8 python3 run195.py cfg/gate_p153c_L128_P0.json gate/res > gate/res1.log 2>&1
  C195_OFFLOAD=1 python3 run195.py cfg/gate_p153c_L128_P0.json gate/res > gate/res2.log 2>&1 ) &
P3=$!
wait $P1 $P2 $P3
echo "gate runs done $(date)"
python3 gate195.py gate/ref/p153c_L128_P0_rand1.json gate/off/p153c_L128_P0_rand1.json gate/res/p153c_L128_P0_rand1.json > logs/gate_$SLURM_JOB_ID.txt 2>&1
tail -3 logs/gate_$SLURM_JOB_ID.txt
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
python3 bench195.py --L 512 --chis 768,1024,1536 --nb 3 --nmv 100 --out logs/bench_L512_$SLURM_JOB_ID.json > logs/bench_L512_$SLURM_JOB_ID.log 2>&1
echo "bench rc=$? $(date)"
kill $SMI

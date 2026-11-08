#!/bin/bash
#SBATCH --job-name=c159
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c159/logs/c159_%j.out
# C159 packed job (Claude Opus agent, 2026-11-07). Resumable stages.
#   A  ring scans of T, K and of the ferro/antiferro neighbourhoods in T and K: 48 tasks (scan/*.jsonl)
#   V  damped VUMPS FES at the C68 anchors of T, K, N: 3 processes (fes/*.json)
#   B  refinement of candidates (ring 16/20, open DMRG 32/64/128): 45 processes (refine/*.json)
cd /scratch/lustre/jkudlerflam/c159
mkdir -p scan refine fes logs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
( for t in $(seq 0 15); do echo "T $t 1500"; done; for t in $(seq 0 15); do echo "K $t 1500"; done
  for t in $(seq 0 5); do echo "TFM $t 1000"; done; for t in $(seq 0 3); do echo "KFM $t 1000"; done
  for t in $(seq 0 5); do echo "TAF $t 1000"; done ) | \
  xargs -P 48 -L 1 sh -c 'timeout 100m python3 scan_tk.py $0 $1 $2 scan > logs/A_$0_$1.log 2>&1'
bash vumps_tk.sh &
timeout 6h python3 refine_tk.py --thr 1.1 --top 8 --procs 45 --chi 128 > logs/B.log 2>&1
wait
python3 summarize_tk.py > summary.txt 2>&1

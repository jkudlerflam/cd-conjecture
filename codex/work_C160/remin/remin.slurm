#!/bin/bash
#SBATCH --job-name=c160re
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c160/remin/logs/remin_%j.out
# Re-minimization of the 48 C160 texture bonds (coordinator, 2026-11-08): N = 2000, PMAX = 300, 8 starts + tiled ring.
cd /scratch/lustre/jkudlerflam/c160/remin
mkdir -p out logs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
( for s in $(seq 0 11); do for l in 0 1 2 3; do echo "$s $l"; done; done ) | \
  xargs -P 48 -L 1 sh -c 'timeout 220m python3 remin.py ../texture/tex_$0.jsonl $1 300 2000 8 out/r_$0_$1.json > logs/r_$0_$1.log 2>&1'
echo done > out/ALLDONE

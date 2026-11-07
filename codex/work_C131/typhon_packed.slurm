#!/bin/bash
#SBATCH --job-name=c131cov
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=16:00:00
#SBATCH --output=logs/c131_%j.out
# Packed version of submit_c131.slurm: all 224 multistart tasks on one node, largest rings first.
cd /scratch/lustre/jkudlerflam/c131_cover
mkdir -p results logs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
( seq 96 127; seq 192 223; seq 64 95; seq 160 191; seq 32 63; seq 128 159; seq 0 31 ) | \
  xargs -P 48 -I{} sh -c 'python3 optimize_q.py --task {} --budget 4000 > logs/task_{}.log 2>&1'
python3 summarize.py results > summary.txt

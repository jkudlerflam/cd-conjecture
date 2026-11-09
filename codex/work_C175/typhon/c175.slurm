#!/bin/bash
#SBATCH --job-name=c175
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c175/logs/c175_%j.out
# C175 packed job (Claude Opus 5.5 agent, 2026-11-09). 14 tasks, 3 BLAS threads each, all at once.
# Effective model (3.2) at phi = pi/2, two-size c at L = 96, 192 (Lb = 32, 64): |w| = 1 below the candidate density (W1),
# |w| = 1.5 and 2 where the inter-valley pair is compact (W15, W2), controls V2 = V3 = 0 (W1M, W2M, MGS limit, c = 2)
# and V3 = 0 (W1R, exact zero-energy resonance, no binding).
cd /scratch/lustre/jkudlerflam/c175
mkdir -p logs
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3
cat typhon/tasks.txt | xargs -P 14 -I CMD sh -c 'CMD >> logs/tasks.log 2>&1'

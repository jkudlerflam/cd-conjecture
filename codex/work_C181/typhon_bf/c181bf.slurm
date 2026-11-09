#!/bin/bash
#SBATCH --job-name=c181bf
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=11:55:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c181/logs_bf/c181bf_%j.out
# C181 backfill job (Claude Opus 5.5 agent, 2026-11-09): the decisive |w| = 1 part of M1 to M6 (phi = 1.53, 1.555,
# n = 3/16, L = 64 and 128 with N +- 1) and the c = 1, c = 2 controls. 16 runs, all at once, 1 thread each, chi to 384.
cd /scratch/lustre/jkudlerflam/c181
mkdir -p logs_bf out_bf
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cat typhon_bf/tasks.txt | xargs -P 16 -I CMD sh -c 'CMD'

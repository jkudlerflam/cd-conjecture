#!/bin/bash
#SBATCH --job-name=c185bf
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c185/bf/logs/c185bf_%j.out
# C185 backfill job (Claude Opus 5.5 agent, 2026-11-09), split from job 17604455 at the coordinator's request.
# 8 blockk runs (kappa = 1 at phi = 1.625, 1.640, 1.655 and the kappa = 0, v = 1 reference, L = 48 and 96), 1 thread each;
# 2 full-chain V/Delta = 8 runs at L = 64 (DM lambda = 0.3 at phi shift 0.06, and the lambda = 0 reference), 2 threads each.
cd /scratch/lustre/jkudlerflam/c185
mkdir -p bf/logs bf/outbk bf/oc
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
(cat typhon/bf_block_tasks.txt | xargs -P 8 -I CMD sh -c 'CMD') &
(export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2; cat typhon/bf_full_tasks.txt | xargs -P 2 -I CMD sh -c 'CMD') &
wait

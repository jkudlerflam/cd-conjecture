#!/bin/bash
#SBATCH --job-name=c186b
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c186b/logs/c186b_%j.out
# C186 backfill job (Claude Opus 5.5 agent, 2026-11-09), split from 17604527 at the coordinator's request.
# Decisive block-DMRG runs of the mixed cusp model (PROOF eq. 3.1): v = 0.3 with the vertex r = sqrt(v), u = 1.035,
# L = 96 and 192 (Lb = 32, 64), mu = onset + 0.03 and + 0.06; plus C173's estimator controls (block_dmrg.py
# unmodified, same configs as C173's T_c1 and T_c2: hard rods c = 1, two decoupled XX legs c = 2) at the same L.
# 8 slots of 2 threads.
cd /scratch/lustre/jkudlerflam/c18x/c186b
mkdir -p logs out
module purge >/dev/null 2>&1; module load anaconda3/2025.06
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2 NUMBA_CACHE_DIR=/scratch/lustre/jkudlerflam/c18x/c186b/numba_cache
cat typhon_bf/tasks.txt | xargs -P 8 -I CMD sh -c 'CMD'
echo done > logs/ALL_DONE

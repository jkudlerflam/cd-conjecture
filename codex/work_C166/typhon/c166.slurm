#!/bin/bash
#SBATCH --job-name=c166
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=320G
#SBATCH --time=40:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c166/logs/c166_%j.out
# C166 packed job (Claude Opus 5.5 agent, 2026-11-08). Two pools on one 48-core node:
#   heavy (4 at a time, 6 BLAS threads each = 24 cores): open chains at L = 512 (4 bonds) and 384 (2 bonds)  -> oc/
#   light (24 at a time, 1 thread each = 24 cores): VUMPS cells n = 2..20 with chi up to 512 (vc/), open chains at
#          L = 128, 256 with observables, pinning-field chains at L = 128, 256 (oc/)
# The job ends with summarize.py > summary.txt.
cd /scratch/lustre/jkudlerflam/c166
mkdir -p logs oc vc
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
( export OPENBLAS_NUM_THREADS=6; cat typhon/tasks_heavy.txt | xargs -P 4 -I CMD sh -c 'CMD >> logs/heavy.log 2>&1' ) &
( export OPENBLAS_NUM_THREADS=1; cat typhon/tasks_light.txt | xargs -P 24 -I CMD sh -c 'CMD >> logs/light.log 2>&1' ) &
wait
python3 summarize.py > summary.txt 2>&1

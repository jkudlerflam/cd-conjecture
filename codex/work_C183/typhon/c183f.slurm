#!/bin/bash
#SBATCH --job-name=c183f
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c183/f/logs/c183f_%j.out
# C183 small job (Claude Opus 5.5 agent, 2026-11-09), the V/Delta = 8 unbound-pair confirmation (coordinator exception
# to the one-job rule).  Points S1 (stratum) and S2 (phi = 1.555) at r = 8, |w| = 1, see typhon/f_points.json.
# Lane A: the two L = 256 open chains, 4 threads each.  Lane B: 4 slots x 2 threads, L = 128 and 64, parity sectors, controls.
cd /scratch/lustre/jkudlerflam/c183/f
mkdir -p logs oc ctrl
( export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4; cat typhon/f_laneA.txt | xargs -P 2 -I CMD sh -c 'CMD' ) &
( export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2; cat typhon/f_laneB.txt | xargs -P 4 -I CMD sh -c 'CMD' ) &
wait

#!/bin/bash
#SBATCH --job-name=c163
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c163/logs/c163_%j.out
# C163 packed job (Claude Opus agent, 2026-11-08). 68 single-core tasks from tasks.txt, 48 at a time:
#   refine_big.py  open-chain DMRG at L = 128, 192, 256 for the C159 points with rising slopes (hi_T_40, dist_T_50,
#                  two 0.03 displacements of hi_T_40) and the slowly falling ones (hi_T_37, hi_K_1)      -> big/
#   fesx.py gen    two-site-cell damped VUMPS FES, chi 8..192, at hi_T_40 and dist_T_50                   -> fes/
#   twisted.py     16 bonds -R_z(phi) + eps X, phi = 1, 2, 2pi/3, pi/2, eps = 0.15, 0.3, DMRG L = 64, 128 -> twisted/
#   path_scan.py   40 attractive-window paths (T, K_I, N), entanglement-peak-first, final DMRG L = 64, 128 -> paths/
cd /scratch/lustre/jkudlerflam/c163
mkdir -p logs big fes twisted paths
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
cat tasks.txt | xargs -P 48 -I CMD sh -c 'CMD >> logs/tasks.log 2>&1'
python3 summarize.py > summary.txt 2>&1

#!/bin/bash
#SBATCH --job-name=c167
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c167/logs/c167_%j.out
# C167 packed job (Claude Opus 5.5 agent, 2026-11-08). 18 tasks, 2 BLAS threads each (36 cores):
# two-site VUMPS along the hi_K_1 harmonic path to chi 192 (256 at lam = 1), the charge-separated bonds
# (q1only, q2only), a four-site cell at lam = 1, lam = 1.01, 1.02 (onset of polarization), and open chains of
# hi_K_1 at L = 128, 256 with edge pinning along the VUMPS order (LRst) plus L = 128 unpinned.
# Digest afterwards locally: python3 xpi_fes.py vc/f_*.json ; python3 ana_path.py vc/*.json
cd /scratch/lustre/jkudlerflam/c167
mkdir -p logs oc vc
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cat typhon/tasks.txt | xargs -P 18 -I CMD sh -c 'CMD >> logs/tasks.log 2>&1'
python3 xpi_fes.py vc/f_*.json > summary.txt 2>&1
python3 ana_path.py vc/*.json >> summary.txt 2>&1

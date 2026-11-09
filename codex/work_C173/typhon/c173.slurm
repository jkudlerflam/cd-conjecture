#!/bin/bash
#SBATCH --job-name=c173
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=40:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c173/logs/c173_%j.out
# C173 packed job (Claude Opus 5.5 agent, 2026-11-09). 12 tasks, 4 BLAS threads each, all at once on one node.
# Effective blockaded model of the Fibonacci cusp (block_dmrg.py) at the symmetric dilute point phi = pi/2,
# mu = -2.2 (n = 0.17) to L = 192, near-symmetric phi = 1.54, 1.52, c = 1 and c = 2 controls at L = 144,
# and the full parity-class qubit chain at V/Delta = 40, 80 (full_dmrg.py).
cd /scratch/lustre/jkudlerflam/c173
mkdir -p logs
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cat typhon/tasks.txt | xargs -P 12 -I CMD sh -c 'CMD >> logs/tasks.log 2>&1'
python3 summ.py > summary.txt 2>&1

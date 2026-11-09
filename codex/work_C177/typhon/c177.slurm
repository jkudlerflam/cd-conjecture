#!/bin/bash
#SBATCH --job-name=c177
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c177/logs/c177_%j.out
# C177 packed job (Claude Opus 5.5 agent, 2026-11-09). 12 slots, 4 BLAS threads each, on one node.
# Full parity-class qubit chain at the Fibonacci cusp, V/Delta = 8, 12, 16, 24 (band-tuned D and h_z, tune.json),
# open chains (fc.py + dmrg2.py) to L = 256 and damped VUMPS (vumps_fc.py) to chi = 384, plus the generic point G8.
cd /scratch/lustre/jkudlerflam/c177
mkdir -p logs oc vc
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cat typhon/tasks.txt | xargs -P 12 -I CMD sh -c 'CMD'
python3 summ177.py . > summary.txt 2>&1

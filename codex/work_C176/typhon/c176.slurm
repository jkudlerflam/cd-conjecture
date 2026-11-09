#!/bin/bash
#SBATCH --job-name=c176
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=320G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c176/logs/c176_%j.out
# C176 packed job (Claude Opus 5.5 agent, 2026-11-09). 38 open-chain runs of the Fibonacci-cusp model (C173 eq. 3.2)
# with the blockade-exact, N-conserving DMRG fibdmrg.py; 24 at a time, 2 BLAS threads each.
cd /scratch/lustre/jkudlerflam/c176
mkdir -p logs out
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cat typhon/tasks.txt | xargs -P 24 -I CMD sh -c 'CMD'
python3 analyze.py out > summary.txt 2>&1

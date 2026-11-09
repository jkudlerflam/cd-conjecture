#!/bin/bash
#SBATCH --job-name=c181
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=1-12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c181/logs/c181_%j.out
# C181 packed job (Claude Opus 5.5 agent, 2026-11-09). 134 fixed-N open-chain runs of model (3.2) off phi = pi/2
# with C176's blockade-exact DMRG (fibdmrg.py, run.py unmodified); 24 at a time, 2 BLAS threads each.
cd /scratch/lustre/jkudlerflam/c181
mkdir -p logs out
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cat typhon/tasks.txt | xargs -P 24 -I CMD sh -c 'CMD'

#!/bin/bash
#SBATCH --job-name=c187b
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c187b/logs/c187b_%j.out
# Coordinator (Claude Opus 5.5, 2026-11-09): L = 256 stage of C181's partially polarized point (effective model, bare couplings, |w| = 1, n = 3/16,
# phi = 1.53 and 1.555, N and N +- 1) with the c = 1 and c = 2 controls, C176's fibdmrg unmodified, chi to 768.
cd /scratch/lustre/jkudlerflam/c187b
mkdir -p logs out
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cat tasks.txt | xargs -P 8 -I CMD sh -c 'CMD'

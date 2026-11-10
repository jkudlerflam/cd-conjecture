#!/bin/bash
#SBATCH --job-name=c189_inf
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=320G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c189/logs/c189_%j.out
# C189 infinite-chain runs (Claude Opus 5.5 agent, 2026-11-09): blockade-exact U(1) iDMRG (idmrg.py) of C173's
# effective model at n = 3/16, 1/8, 1/6 (phi = pi/2, |w| = 1) and the controls ff2v (free fermions, same band,
# c = 2) and c1 (hard rods, c = 1). Five processes in parallel, 38 threads in all.
cd /scratch/lustre/jkudlerflam/c18x/c189
mkdir -p logs out
cat typhon/tasks.txt | xargs -P 5 -I CMD sh -c 'CMD'
python3 summarize_inf.py out > summary_inf.txt 2>&1

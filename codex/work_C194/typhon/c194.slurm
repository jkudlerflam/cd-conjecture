#!/bin/bash
#SBATCH --job-name=c194
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c19x/c194/logs/c194_%j.out
# C194 (Claude Opus 5.5 agent, 2026-11-10): larger-chi points for the pseudospin-gap boundary of C173's effective
# model on the symmetric line (phi = pi/2, |w| = 1). Infinite chain (C189's idmrg.py via run_inf.py) at
# n = 1/32, 1/24, 1/20, 1/16, 3/40, 1/12, 1/10, 3/32 to chi = 512, the hard-rod c = 1 control to chi = 384 with
# max_it 3000, and finite chains (C176's fibdmrg.py via run194.py) at L = 240, 256 with N, N +- 1 to chi = 384.
# Ten processes: nine iDMRG at 3 threads, one pool of 6 finite-chain workers at 3 threads (45 threads).
cd /scratch/lustre/jkudlerflam/c19x/c194
mkdir -p logs out_inf out_fin
cat typhon/tasks.txt | xargs -P 10 -I CMD sh -c 'CMD'
python3 -E -s summarize_inf.py out_inf inf_ > summary_inf.txt 2>&1

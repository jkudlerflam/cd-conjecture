#!/bin/bash
#SBATCH --job-name=c198
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=360G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c19x/c198/logs/c198_%j.out
# C198 (coordinator, 2026-11-10): the size sequence of open-chain c_two in C173's effective model at n = 3/16 on the
# symmetric line, L = 128 to 1024, with C176's code unchanged (run.py, fibdmrg.py, analyze.py). L = 128, 192, 256
# resume from C176's checkpoints (copied, read only at source) and add chi = 1536. Controls c1 (hard rods) and c2
# (two XX legs) at L = 512 and 1024. Nine runs at once, 44 BLAS threads in total. A resubmission resumes every run.
cd /scratch/lustre/jkudlerflam/c19x/c198
mkdir -p logs out
cat typhon/tasks.txt | xargs -P 9 -I CMD sh -c 'CMD'
python3 analyze.py out > summary.txt 2>&1

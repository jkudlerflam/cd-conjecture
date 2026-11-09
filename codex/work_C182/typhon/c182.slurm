#!/bin/bash
#SBATCH --job-name=c182
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=320G
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c182/logs/c182_%j.out
# C182 packed job (Claude Opus 5.5 agent, 2026-11-09). Stage 1: cluster ladders (clusters.py, ladder6.py) and the
# DMRG phase-diagram screen of the Fibonacci-cusp model (C176's fibdmrg.py, unmodified). Stage 2: converged series
# at the points chosen by select_stage2.py (rule registered in RESULT.md). 24 tasks at a time, 2 BLAS threads each.
cd /scratch/lustre/jkudlerflam/c182
mkdir -p logs out out2
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cat typhon/tasks1.txt | xargs -P 24 -I CMD sh -c 'CMD'
python3 select_stage2.py > logs/select_stage2.log 2>&1
cat typhon/tasks2.txt | xargs -P 24 -I CMD sh -c 'CMD'
python3 analyze.py out > summary1.txt 2>&1
python3 analyze.py out2 > summary2.txt 2>&1

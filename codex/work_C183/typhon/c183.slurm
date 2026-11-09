#!/bin/bash
#SBATCH --job-name=c183
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c183/logs/c183_%j.out
# C183 packed job (Claude Opus 5.5 agent, 2026-11-09). One node, 48 cores.
# Phases 0-3: 24 workers x 2 BLAS threads (tuning, density ladders, L = 32/64 screen, free controls).
# Phase 4: 16 workers x 3 threads (L = 128/256 open chains, parity sectors, VUMPS FES, controls) on the points
# selected by the rule registered in RESULT.md.
cd /scratch/lustre/jkudlerflam/c183
mkdir -p logs
python3 pipeline183.py --workers 24 --threads 2 --phase4-workers 16 --phase4-threads 3

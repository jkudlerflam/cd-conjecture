#!/bin/bash
#SBATCH --job-name=c187d
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c187d/logs/c187d_%j.out
# Coordinator (Claude Fable 5, 2026-11-09): continuation of the C187b L = 256 runs and their L = 128 partners to chi = 1024.
# Fallback for the overlap steps. A tag that is complete or running elsewhere is skipped.
cd /scratch/lustre/jkudlerflam/c187d
bash body.sh tasks_all.txt 10

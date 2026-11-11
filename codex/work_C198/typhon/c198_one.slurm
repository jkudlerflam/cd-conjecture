#!/bin/bash
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c19x/c198/logs/%x_%j.out
# C198, one run per job so that each can start in a gap of the queue (coordinator, 2026-11-10).
# Usage: sbatch --job-name=c198_<tag> --cpus-per-task=<t> --mem=<m> typhon/c198_one.slurm <tag>
cd /scratch/lustre/jkudlerflam/c19x/c198
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
timeout 4300m python3 run.py typhon/cfg/$1.json out > logs/$1.log 2>&1
echo "rc=$? $(date)"

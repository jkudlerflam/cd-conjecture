#!/bin/bash
#SBATCH --job-name=c143_srch
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c143/slurm_%j.out
cd /scratch/lustre/jkudlerflam/c143
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/opt/sns/anaconda3/2024.06/bin/python3 -B typhon_search.py /scratch/lustre/jkudlerflam/c143/search_v2_10_14.jsonl 48 600

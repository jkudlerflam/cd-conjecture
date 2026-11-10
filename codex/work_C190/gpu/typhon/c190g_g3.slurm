#!/bin/bash
#SBATCH --job-name=c190g_g3
#SBATCH --partition=apollo
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c190g/logs/c190g_g3_%j.out
# C190g gate G3 (Claude Opus 5.5 agent, 2026-11-09): the L = 64 local pilot (collaborator's field, pilot schedule)
# on the GPU.  Usage: sbatch c190g_g3.slurm <P>
cd /scratch/lustre/jkudlerflam/c18x/c190g
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export C190G_BACKEND=cupy C190G_PROFILE=1
mkdir -p g3out
python3 run190g.py cfg_g3/p153c_L64_P$1.json g3out > logs/G3_P$1.log 2>&1

#!/bin/bash
#SBATCH --job-name=c168b
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c168/b/logs/c168b_%j.out
# C168 follow-up: FES (two-site and four-site cells, chi to 256) and DMRG reseeds (L = 128, two seeds; L = 192)
# for the two stage-A points that met the first stage of the pre-fixed rule (10193, 450) and one control (10540).
cd /scratch/lustre/jkudlerflam/c168/b
mkdir -p logs fes dm
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
xargs -P 15 -I{} sh -c '{}' < tasks.txt
python3 fesfit_c158.py fes/*.json > fesfit.out 2>&1

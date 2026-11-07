#!/bin/bash
#SBATCH -J c137win
#SBATCH -p typhon
#SBATCH -N 1
#SBATCH -c 48
#SBATCH --mem=180G
#SBATCH -t 24:00:00
#SBATCH -o /scratch/lustre/jkudlerflam/c137_window/slurm-%j.out
# C137: random-bond window survey n = 3..13, rings to L = 16, main bonds to n = 18. One node, packed with xargs.
cd /scratch/lustre/jkudlerflam/c137_window/work_C137
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 C137_SPLIT=1
mkdir -p results/rings logs
python3 tasks_typhon.py > tasks.txt
cat tasks.txt | xargs -P 46 -L 1 sh -c 'python3 $0 $1 $2 $3 > logs/$0_$1_$2.log 2>&1'

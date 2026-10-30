#!/bin/bash
# open-chain block-sparse DMRG profiles, packed into one allocation (per-user job cap)
#SBATCH -J n2prof
#SBATCH -p typhon
#SBATCH -c 48
#SBATCH --mem=180G
#SBATCH -t 3-00:00:00
#SBATCH -o logs/%x_%j.out
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd $SLURM_SUBMIT_DIR
case ${1:-A} in
A) T="M0 32 1024|M1 32 1024|M2 32 1024|C0 32 512|C1 32 512|C0 48 512|C1 48 512";;
B) T="M0 48 1024|M1 48 1024|M2 48 1024|C0 64 512|C1 64 512|C0 32 256|C1 32 256";;
C) T="M0 32 1536|M1 32 1536|M2 32 1536|M0 32 512|M1 32 512|M2 32 512";;
esac
echo "$T" | tr '|' '\n' | xargs -P 6 -n 3 sh -c 'python3 -B run_profile.py $0 $1 $2 > logs/prof_$0_L$1_chi$2.log 2>&1'

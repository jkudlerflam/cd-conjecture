#!/bin/bash
# VUMPS finite-entanglement runs, packed into one allocation (per-user job cap).
# Stage 1 builds the payload-only initial states, stage 2 runs the models.
#SBATCH -J vC44n2fes
#SBATCH -p typhon
#SBATCH -c 48
#SBATCH --mem=180G
#SBATCH -t 3-00:00:00
#SBATCH -o logs/%x_%j.out
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd $SLURM_SUBMIT_DIR
mkdir -p results/init_cache
for cp in 32 40 48 64 80 96 128 160; do echo $cp; done |
  xargs -P 8 -n 1 sh -c 'python3 -B -c "import run_fes; run_fes.payload_state($0)" > logs/cache_$0.log 2>&1'
{ for chi in 288 240 192 144 120 96; do for m in M0 M1 M2; do echo "$m $chi"; done; done
  for chi in 160 128 96 64 48 32; do for m in C0 C1; do echo "$m $chi"; done; done; } |
  xargs -P 12 -n 2 sh -c 'python3 -B run_fes.py $0 $1 3 > logs/fes_$0_$1.log 2>&1'

#!/bin/bash
#SBATCH --job-name=c186
#SBATCH --partition=typhon
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/lustre/jkudlerflam/c18x/c186/logs/c186_%j.out
# C186 packed job (Claude Opus 5.5 agent, 2026-11-09). 24 slots of 2 threads.
# (A) pairbind2.py: binding of two dressed quasiparticles of the mixed cusp model (PROOF eq. 3.1) on L = 48 rings,
#     N_max = 4 and 5, u = 0.522, 1.035, 1.531, 2 (valley momenta on the grid), v / v_c = 0.25 ... 2, r = sqrt(v);
#     r = 0 calibrations at N_max = 2.
# (B) bm_driver.py: open-chain block DMRG of the mixed model at u = 1.035, L = 96 and 192, two densities, with
#     r = 0 controls (v = 0 bound pair, v = 0.4 unbound).
cd /scratch/lustre/jkudlerflam/c18x/c186
mkdir -p logs out
module purge >/dev/null 2>&1; module load anaconda3/2025.06
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2 NUMBA_CACHE_DIR=/scratch/lustre/jkudlerflam/c18x/c186/numba_cache
cat typhon/tasks.txt | xargs -P 24 -I CMD sh -c 'CMD'
echo done > logs/ALL_DONE

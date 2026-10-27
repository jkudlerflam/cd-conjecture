#!/bin/bash
#SBATCH --job-name=d2fu_obcP
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=32
#SBATCH --mem=80G
#SBATCH --time=3-00:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/obcP_%A_%a.out
# Open-chain DMRG (fdmrg.py): array task i runs lines 4i+1..4i+4 of obc_tasks.txt in parallel,
# 8 BLAS threads each.  Line format: model[:edge] L seed chis...   (":edge" = full field on the end
# sites, the convention of d2_classification/dmrg.py; default = half field, sum of bond terms).
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
for k in 1 2 3 4; do
  LINE=$((4*SLURM_ARRAY_TASK_ID+k))
  read MODEL L SEED CHIS <<< "$(sed -n "${LINE}p" obc_tasks.txt)"
  EXTRA=""; PRE="obc"
  if [[ "$MODEL" == *:edge ]]; then MODEL=${MODEL%:edge}; EXTRA="--edgefield"; PRE="obcedge"; fi
  python3 fdmrg.py "$MODEL" --L "$L" --seed "$SEED" --chis $CHIS --sweeps 4 --etol 1e-10 --win 64 --rmax 80 $EXTRA \
      --out "data/${PRE}_${MODEL}_L${L}_s${SEED}.json" > "logs/${PRE}_${MODEL}_L${L}_s${SEED}.log" 2>&1 &
done
wait

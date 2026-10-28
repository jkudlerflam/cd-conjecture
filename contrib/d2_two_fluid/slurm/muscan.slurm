#!/bin/bash
#SBATCH --job-name=tf_mu
#SBATCH --partition=typhon
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/mu_%A_%a.out
# iDMRG chemical-potential scan.  usage: sbatch --array=0-N muscan.slurm NAME 'p:t,phi,V,{MU},gL,|gR|,arg' MU0 DMU CHIS...
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
NAME=$1; TPL=$2; MU0=$3; DMU=$4; shift 4
MU=$(python3 -c "print(round($MU0+$DMU*$SLURM_ARRAY_TASK_ID,6))")
SPEC=${TPL/\{MU\}/$MU}
python3 run_idmrg.py "$SPEC" --chis "$@" --steps 300 --tol 1e-9 --out data/mu_${NAME}.jsonl --label "$MU"

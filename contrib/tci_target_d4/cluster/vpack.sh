#!/bin/bash
#SBATCH --partition=typhon
#SBATCH --qos=short
#SBATCH --output=logs/%x_%j.out
# usage: sbatch --job-name=X --cpus-per-task=C --mem=M --time=T vpack.sh TASKFILE FIRST LAST THREADS
# runs lines FIRST..LAST of TASKFILE (arguments of vumps_run.py) concurrently
cd /scratch/lustre/jkudlerflam/tci_d4/tci_target_d4
TF=$1; FIRST=$2; LAST=$3; TH=$4
export OPENBLAS_NUM_THREADS=$TH OMP_NUM_THREADS=$TH MKL_NUM_THREADS=$TH
for i in $(seq $FIRST $LAST); do
  LINE=$(sed -n "$((i+1))p" "$TF")
  python3 -u vumps_run.py $LINE > logs/$(basename $TF .txt)_$i.log 2>&1 &
done
wait

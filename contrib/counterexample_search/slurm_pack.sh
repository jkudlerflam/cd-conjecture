#!/bin/bash
#SBATCH --job-name=cdpack
#SBATCH --partition=typhon
#SBATCH --output=logs/%x_%j.out
# usage: sbatch --cpus-per-task=C --mem=M --time=T slurm_pack.sh SCRIPT TASKFILE FIRST LAST THREADS
# runs lines FIRST..LAST of TASKFILE concurrently, each with THREADS BLAS threads
cd /scratch/lustre/jkudlerflam/cdsearch
SCRIPT=$1; TF=$2; FIRST=$3; LAST=$4; TH=$5
export OPENBLAS_NUM_THREADS=$TH OMP_NUM_THREADS=$TH MKL_NUM_THREADS=$TH
for i in $(seq $FIRST $LAST); do
  LINE=$(sed -n "$((i+1))p" "$TF")
  echo "start $i: $LINE"
  python3 -u $SCRIPT $LINE > logs/pack_$(basename $TF .txt)_$i.log 2>&1 &
done
wait
echo all done

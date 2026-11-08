#!/bin/bash
export PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /scratch/lustre/jkudlerflam/c153
LINE="$1"
NAME=$(echo "$LINE" | tr ' /' '__' | cut -c1-150)
python3 $LINE > logs/$NAME.log 2>&1
echo "finished $LINE"

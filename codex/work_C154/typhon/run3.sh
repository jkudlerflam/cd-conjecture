#!/bin/bash
cd /scratch/lustre/jkudlerflam/c154
PY=/opt/sns/anaconda3/2024.06/bin/python3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
set -- $1
kind=$1; shift
if [ "$kind" = "grid" ]; then
  $PY -s -B cyc_grid.py "$@" >> cyc_grid.out 2>&1
else
  $PY -s -B path_dmrg.py 1282 "$1" "$2" "$3" "$4" >> path1282_L$1.out 2>&1
fi

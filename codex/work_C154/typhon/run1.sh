#!/bin/bash
cd /scratch/lustre/jkudlerflam/c154
set -- $1
PY=/opt/sns/anaconda3/2024.06/bin/python3
if [ "$1" = "ed" ]; then
  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 $PY -s -B ed_sec.py "$2" "$3" "$4" "$5" "$6" > "ed24_${3//,/_}.out" 2>&1
else
  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 $PY -s -B dmrg_cyc.py "$2" "$3" "$4" "$5" > "${4%.jsonl}.out" 2>&1
fi

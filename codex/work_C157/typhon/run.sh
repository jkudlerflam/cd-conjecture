#!/bin/bash
cd /scratch/lustre/jkudlerflam/c157
set -- $1
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 /opt/sns/anaconda3/2024.06/bin/python3 -s -B dmrg_run.py "$1" "$2" "$3" "$4" > "${3%.jsonl}.out" 2>&1

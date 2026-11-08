#!/bin/bash
# args: L chis label t phi Delta V mu
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /scratch/lustre/jkudlerflam/c152/j2
python3 dmrg_pt.py $1 $2 dmrg_$1.jsonl $3 $4 $5 $6 $7 $8 > /dev/null 2>> dmrg_err.log

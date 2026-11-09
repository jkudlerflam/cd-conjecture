#!/bin/bash
# one task line: KIND L chis label args...   KIND = pt (parity-class point, spin couplings) or ctrl
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /scratch/lustre/jkudlerflam/c18x/c184
kind=$1; L=$2; chis=$3; lab=$4; shift 4
if [ "$kind" = "pt" ]; then
  python3 dmrg_pt.py $L $chis out/dmrg_$L.jsonl $lab "$@" > /dev/null 2>> out/err.log
else
  python3 controls.py $L $chis out/ctrl_$L.jsonl $lab > /dev/null 2>> out/err.log
fi

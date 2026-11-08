#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
cd "codex/work_C152"
python3 dmrg_pt.py $1 $2 local_dmrg/dmrg_$1.jsonl $3 $4 $5 $6 $7 $8 >> local_dmrg/log.txt 2>&1

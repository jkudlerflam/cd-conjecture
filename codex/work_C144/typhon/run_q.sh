#!/bin/bash
# one Q-scan task: $1 = task line, output appended per label
cd /scratch/lustre/jkudlerflam/c144
lab=$(echo "$1" | cut -d: -f1)
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /opt/sns/anaconda3/2024.06/bin/python3 -s -B qscan.py 10 "$1" > q_$lab.jsonl 2> q_$lab.err

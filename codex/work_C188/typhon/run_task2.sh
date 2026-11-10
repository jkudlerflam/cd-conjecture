#!/bin/bash
# C188 job 2 dispatcher. ring30 pops the next point from ring30_queue.txt (300 Lanczos steps, sectors 15 = S, 10 = F, 14 = current).
cd /scratch/lustre/jkudlerflam/c18x/c188
P="/opt/sns/anaconda3/2024.06/bin/python3 -s -B"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
set -- $1
k=$1; shift
case $k in
ring24) OMP_NUM_THREADS=1 $P ringk_run.py 24 "$1" $(seq -s, 0 12) 200 out2/ring24.jsonl ./ringk >> out2/ring24.log 2>&1 ;;
ring30) p=$(head -1 ring30_queue.txt); [ -z "$p" ] && exit 0; sed -i 1d ring30_queue.txt; OMP_NUM_THREADS=8 $P ringk_run.py 30 "$p" 15,10,14 300 out2/ring30.jsonl ./ringk >> out2/ring30.log 2>&1 ;;
dmrg) OMP_NUM_THREADS=1 $P dmrg_tau.py "$1" "$2" out2/dmrg_L$1.jsonl "$3" >> out2/dmrg.log 2>&1 ;;
esac

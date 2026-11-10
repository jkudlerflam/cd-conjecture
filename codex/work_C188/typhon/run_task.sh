#!/bin/bash
cd /scratch/lustre/jkudlerflam/c18x/c188
P="/opt/sns/anaconda3/2024.06/bin/python3 -s -B"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
set -- $1
k=$1; shift
case $k in
ring24) OMP_NUM_THREADS=1 $P ringk_run.py 24 "$1" $(seq -s, 0 12) 200 out/ring24.jsonl ./ringk >> out/ring24.log 2>&1 ;;
ring30) p=$(head -1 ring30_queue.txt); [ -z "$p" ] && exit 0; sed -i 1d ring30_queue.txt; OMP_NUM_THREADS=8 $P ringk_run.py 30 "$p" 15,10 300 out/ring30_v2.jsonl ./ringk >> out/ring30.log 2>&1 ;;
dmrg) OMP_NUM_THREADS=1 $P dmrg_tau.py "$1" "$2" out/dmrg_L$1.jsonl "$3" >> out/dmrg.log 2>&1 ;;
esac
#######################
exit 0
# C188 task dispatcher, revised 2026-11-09 ~14:50 while the job ran. ring30 tasks ignore their argument and pop the
# next point from ring30_queue.txt, 300 Lanczos steps in sectors 15 (S) and 10 (F), because 130 steps did not converge
# at L = 30. The 'exit 0' sits at byte 626, the length of the original script, so that bash instances already running
# the old script (which resume reading at that offset when their command returns) exit cleanly.

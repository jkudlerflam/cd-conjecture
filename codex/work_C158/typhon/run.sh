#!/bin/bash
# usage: run.sh "<kind> <tag> <command...>"
cd /scratch/lustre/jkudlerflam/c158
line="$1"; kind=${line%% *}; rest=${line#* }; tag=${rest%% *}; cmd=${rest#* }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PATH=/opt/sns/anaconda3/2024.06/bin:$PATH
echo "start $tag $(date)" >> logs/status.txt
bash -c "$cmd" > logs/$tag.log 2>&1
echo "end $tag rc=$? $(date)" >> logs/status.txt

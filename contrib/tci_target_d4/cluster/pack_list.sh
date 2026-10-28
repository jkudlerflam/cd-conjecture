#!/bin/bash
#SBATCH --partition=typhon
#SBATCH --qos=short
#SBATCH --output=logs/%x_%j.out
# usage: sbatch ... pack_list.sh NCHUNK_TOTAL POINTS OUTPREFIX "c1 c2 ..." SIZES...   (reruns listed chunks)
cd /scratch/lustre/jkudlerflam/tci_d4/tci_target_d4
NT=$1; PTS=$2; OUT=$3; LIST=$4; shift 4
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
for i in $LIST; do
  python3 -u chart.py --sizes "$@" --points $PTS --chunk $i $NT --nlev 8 --out data/${OUT}r_c$i.jsonl > logs/${OUT}r_c$i.log 2>&1 &
done
wait

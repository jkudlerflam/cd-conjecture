#!/bin/bash
#SBATCH --partition=typhon
#SBATCH --qos=short
#SBATCH --output=logs/%x_%j.out
# usage: sbatch --job-name=X --cpus-per-task=C --mem=M --time=T pack.sh FIRSTCHUNK NCHUNK_TOTAL NPROC POINTS OUTPREFIX SIZES...
# runs NPROC chunks (FIRSTCHUNK .. FIRSTCHUNK+NPROC-1 of NCHUNK_TOTAL) of chart.py concurrently, 1 BLAS thread each
cd /scratch/lustre/jkudlerflam/tci_d4/tci_target_d4
FIRST=$1; NT=$2; NP=$3; PTS=$4; OUT=$5; shift 5
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
for i in $(seq $FIRST $((FIRST+NP-1))); do
  python3 -u chart.py --sizes "$@" --points $PTS --chunk $i $NT --nlev 8 --out data/${OUT}_c$i.jsonl > logs/${OUT}_c$i.log 2>&1 &
done
wait
echo all done

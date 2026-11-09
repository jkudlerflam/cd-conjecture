#!/bin/bash
# Quick test of the Fibonacci-cusp effective model at the V/Delta = 8 fitted couplings (see README.txt).
# Usage: bash run_desktop.sh [parallel_jobs]   (default: number of cores / 2)
cd "$(dirname "$0")"
mkdir -p out logs
NC=$( (sysctl -n hw.ncpu 2>/dev/null || nproc) )
P=${1:-$(( NC / 2 ))}; [ "$P" -lt 1 ] && P=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2
echo "running 32 tasks, $P at a time; smallest sizes first. Progress: tail -f logs/*.log ; results: python3 summarize.py"
cat tasks.txt | xargs -P "$P" -I CMD sh -c 'CMD'
python3 summarize.py

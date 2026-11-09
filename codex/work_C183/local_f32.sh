#!/bin/bash
# C183: stage F points S1, S2 at L = 32 locally (same chis as the Typhon L = 64 runs), both parity sectors.
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p local/f
for t in S1 S2; do
  python3 fc.py $t 32 --cp @typhon/$t.json --chis 32 64 96 128 192 192 192 --outdir local/f | tail -1
  for p in odd even; do
    python3 fc.py ${t}$p 32 --cp @typhon/$t.json --chis 32 64 96 128 160 160 --outdir local/f --parity $p --nocorr | tail -1
  done
done

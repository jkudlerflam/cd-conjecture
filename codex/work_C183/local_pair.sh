#!/bin/bash
# C183 local pilot: unbound (r = 5) against bound (r = 12) at |w| = 1, phi = pi/2, n_d ~ 0.12.
# Open chains L = 32, 64 and the parity sectors at the same L.
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
for spec in "5 0.45 r5" "12 0.04 r12"; do
  set -- $spec
  CP=$(python3 -c "import json,tune183 as TU; T=json.load(open('local/p$1_tune.json')) if False else TU.tune($1,1.0,1.5707963267948966); print(json.dumps(TU.couplings(T,$2)))")
  echo "$3 $CP"
  for L in 32 64; do
    python3 fc.py P$3 $L --cp "$CP" --chis 32 64 96 128 160 160 --outdir local | tail -1
    for par in odd even; do
      python3 fc.py P${3}_$par $L --cp "$CP" --chis 32 64 96 128 160 160 --outdir local --parity $par --nocorr | tail -1
    done
  done
done

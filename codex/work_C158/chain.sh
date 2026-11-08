#!/bin/bash
# C158 local chaining: run fesx.py one bond dimension at a time (each call capped at 290 s), restarting from the npz.
# usage: chain.sh TAG "MODEL ARGS" "CHI1 CHI2 ..."
cd "$(dirname "$0")"
tag=$1; margs=$2; chis=$3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
for c in $chis; do
  timeout 290 python3 fesx.py $margs --chis $c --maxit 800 --out local/$tag.json --npz local/$tag.npz >> local/$tag.log 2>&1
  rc=$?; echo "chi=$c rc=$rc $(date +%T)" >> local/$tag.log
  [ $rc -ne 0 ] && break
done

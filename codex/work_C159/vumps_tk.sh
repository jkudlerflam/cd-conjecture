#!/bin/bash
# C159 stage V: FES (damped VUMPS, two-site cell) at the C68 anchors of T, K and N (method controls in the
# strata; a one- or three-site cell cannot hold an injective Theta-symmetric MPS in T, by LSM, and does not converge).
mkdir -p fes logs
run() { tag=$1; shift; timeout 5h python3 fesx.py gen "$@" --out fes/$tag.json > logs/V_$tag.log 2>&1; }
T1="1.1 0.5 -0.2 -0.5 0.9 0.1 0.2 -0.1 0.4"
run T_C68 --J $T1 --B 0 0 0 --n 2 --chis 8 12 16 24 32 48 64 96 128 &
run K_C68 --J 1.1 0.5 0 -0.5 0.9 0 0 0 0.4 --B 0.1 0.2 0 --n 2 --chis 8 12 16 24 32 48 64 96 128 &
run N_C68 --J $T1 --B 0.3 0.4 0.5 --n 2 --chis 8 12 16 24 32 48 64 96 128 &
wait

#!/bin/bash
# supersymmetric t-J as a bosonic chain: h = -P_g on C^{1|2} (hole bosonic, electrons fermionic),
# canonical sectors [N_hole, N_up] with N_down = N_e - N_up; chemical potential added in analysis.
for L in 4 6 8 10 12 14; do
  for Ne in $(seq 0 $L); do
    Nup=$(( (Ne + 1) / 2 )); Nh=$(( L - Ne ))
    python3 run.py --model graded --par -1 0 1 1 0 0 0 --L $L --nlev 3 --targets "[[$Nh, $Nup]]" --out data/tj.jsonl > /dev/null
  done
done
for Ne in 8 12; do Nup=$(( Ne / 2 )); python3 run.py --model graded --par -1 0 1 1 0 0 0 --L 16 --nlev 3 --ks near --targets "[[$(( 16 - Ne )), $Nup]]" --out data/tj.jsonl > /dev/null; done

#!/bin/bash
# d = 3 family scans that fit on a laptop (L <= 12).
cd "$(dirname "$0")"
O=data/families_d3.jsonl
for n in 0.3333 0.5 0.6667 0.8333; do python3 canonical.py --family tj --dens $n --Ls 6 12 --out $O; done
python3 canonical.py --family tj --dens 0.75 --Ls 8 12 --out $O
python3 canonical.py --family tj --dens 0.5 --Ls 4 8 12 --out $O
python3 canonical.py --family tj --dens 0.25 --Ls 4 8 12 --out $O
for s in 1; do for n in 0.5 0.6667; do python3 canonical.py --family tj --par sign=$s --dens $n --Ls 6 12 --out $O; done; done
for th in 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8 4.0 4.2 4.4 4.6 4.8 5.0 5.2 5.4 5.6 5.8 6.0; do
  for n in 0.3333 0.6667; do python3 canonical.py --family osp12 --par theta=$th --dens $n --Ls 6 9 12 --out $O; done
done
for th in 0.6 0.785398 0.9 1.0 1.1 1.2 1.3 1.4 1.5; do python3 scan.py --family blbq --par theta=$th --Ls 6 9 12 --out $O; done
for D in -0.5 0.0 0.5 1.0; do for dl in 0.0 0.5 1.0; do python3 scan.py --family spin1xxzd --par delta=$dl D=$D --Ls 6 8 10 12 --out $O; done; done

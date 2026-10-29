#!/bin/bash
# d = 4 family scans at L <= 8 (full charge scan) for a first look.
cd "$(dirname "$0")"
O=data/families_d4_small.jsonl
for x in 0.25 0.5 0.75 1.0 1.5; do python3 scan.py --family kk --par x=$x y=$x --Ls 4 8 --out $O; done
for x in 0.0 0.25 0.5; do for y in 0.75 1.0; do python3 scan.py --family kk --par x=$x y=$y --Ls 4 8 --out $O; done; done
for lam in 0.0 0.5 1.0 2.0 4.0 8.0; do python3 scan.py --family twochain --par lam=$lam mode=zz --Ls 4 8 --out $O; done
for lam in 0.25 0.5 1.0 -0.5 -1.0; do python3 scan.py --family twochain --par lam=$lam mode=cross --Ls 4 8 --out $O; done
for lam in 0.5 1.0 2.0 -0.5 -1.0 -2.0; do python3 scan.py --family twochain --par lam=$lam mode=rung --Ls 4 8 --out $O; done
for U in 0 1 2 4 8 16; do for n in 0.5 1.0; do python3 canonical.py --family hubbard --par U=$U --dens $n --Ls 4 8 --out $O; done; done
for U in -1 -2 -4; do python3 canonical.py --family hubbard --par U=$U --dens 1.0 --Ls 4 8 --out $O; done
# graded gl(3|1): bosons 0,1,2 antiferro (sign +1), fermion 3; densities over species 1,2,3
for nf in 0.25 0.5; do nb=$(python3 -c "print((1-$nf)/3)"); python3 canonical.py --family gperm --par "par=[0,0,0,1]" sign=1 --dens $nb $nb $nf --Ls 4 8 --out $O; done
# gl(1|3): one boson, three fermions, sign -1 makes the fermion SU(3) part antiferromagnetic
for nb in 0.25 0.5; do nf=$(python3 -c "print((1-$nb)/3)"); python3 canonical.py --family gperm --par "par=[0,1,1,1]" sign=-1 --dens $nf $nf $nf --Ls 4 8 --out $O; done
# gl(2|2): both signs, half bosons
python3 canonical.py --family gperm --par "par=[0,0,1,1]" sign=1 --dens 0.25 0.25 0.25 --Ls 4 8 --out $O
python3 canonical.py --family gperm --par "par=[0,0,1,1]" sign=-1 --dens 0.25 0.25 0.25 --Ls 4 8 --out $O
python3 canonical.py --family gperm --par "par=[0,0,1,1]" sign=1 --dens 0.375 0.125 0.125 --Ls 8 --out $O

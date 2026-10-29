#!/bin/bash
cd "$(dirname "$0")"
O=data/families_d4_L12.jsonl
python3 scan.py --family kk --par x=0.25 y=0.25 --Ls 4 8 12 --nlev 2 --out $O
python3 canonical.py --family gperm --par 'par=[0,0,0,1]' sign=1 --dens 0.25 0.25 0.25 --Ls 4 8 12 --out $O
python3 scan.py --family kk --par x=0.5 y=0.5 --Ls 4 8 12 --nlev 2 --out $O
python3 canonical.py --family hubbard --par U=4 --dens 1.0 --Ls 4 8 12 --out $O
python3 scan.py --family twochain --par lam=0.5 mode=cross --Ls 4 8 12 --nlev 2 --out $O
python3 canonical.py --family hubbard --par U=4 --dens 0.5 --Ls 4 8 12 --out $O
python3 scan.py --family kk --par x=0.0 y=1.0 --Ls 4 8 12 --nlev 2 --out $O

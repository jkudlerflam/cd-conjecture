#!/usr/bin/env python3
"""Pick the points followed with FES and open chains, from data/ring_summary.json.
 (i)  every stable flag (extra neutral level at q = 0 with x0 < 1.3, steady in L), lowest x0 first;
 (ii) per cut, liquid points that touch an Ising-like point (eight-neighbour), the junction neighbourhoods,
      ranked by cS(20), at most --nj per cut;
 (iii) the C68 anchor and point 1282 as c = 1 references.
Writes tasks_fes.txt and tasks_obc.txt.
Written for contrib/d2_floating_ising (Claude Opus 5.5 subagent, 2026-10-29).
"""
import os, sys, json, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import models_fi as MF
ap = argparse.ArgumentParser()
ap.add_argument('--nflag', type=int, default=12)
ap.add_argument('--nj', type=int, default=3)
ap.add_argument('--chis', default='8 12 16 24 32 48 64 96')
a = ap.parse_args()
S = json.load(open(os.path.join(HERE, 'data', 'ring_summary.json')))


def name(k):
    s = S[k]
    return s['fam'] + ':' + ','.join(f'{x}={v:.6g}' for x, v in sorted(s['p'].items()))


pick = []
fl = sorted([k for k in S if S[k]['stable']], key=lambda k: S[k]['L20']['x0'])
pick += [(k, 'flag') for k in fl[:a.nflag]]
for cut in MF.CUTS:
    cand = []
    for k, s in S.items():
        if not k.startswith(cut + '_') or not s['gapless'] or s['L20']['cS'] > 1.5:
            continue
        i, j = int(k.split('_')[1]), int(k.split('_')[2])
        nb = [f'{cut}_{i+di:02d}_{j+dj:02d}' for di in (-1, 0, 1) for dj in (-1, 0, 1) if di or dj]
        if any(n in S and S[n]['ising'] for n in nb):
            cand.append(k)
    cand.sort(key=lambda k: -S[k]['L20']['cS'])
    pick += [(k, 'junction') for k in cand[:a.nj] if k not in dict(pick)]
fes, obc = [], []
for k, why in pick:
    nm = name(k)
    fes.append(f"OPENBLAS_NUM_THREADS=8 python3 fes_fi.py '{nm}' --tag {k} --cells 1 2 3 4 5 6 7 8 9 10 11 12 --chis {a.chis}")
    obc.append(f"OPENBLAS_NUM_THREADS=8 python3 obc_fi.py '{nm}' --tag {k} --L 64 --chis 32 64 128")
    obc.append(f"OPENBLAS_NUM_THREADS=8 python3 obc_fi.py '{nm}' --tag {k} --L 128 --chis 32 64 128")
    print(k, why, nm, 'x0(16,18,20)=', [round(S[k][f'L{L}']['x0'], 3) for L in (16, 18, 20)], 'cS20=', round(S[k]['L20']['cS'], 3))
for nm, tg in (('Z2:g=0,u=0.4,D=0.5,b=0', 'anchor'), ('P1282:s=0,t=0', 'p1282')):
    fes.append(f"OPENBLAS_NUM_THREADS=8 python3 fes_fi.py '{nm}' --tag {tg} --cells 1 2 3 4 5 6 7 8 9 10 11 12 --chis {a.chis}")
open(os.path.join(HERE, 'tasks_fes.txt'), 'w').write('\n'.join(fes) + '\n')
open(os.path.join(HERE, 'tasks_obc.txt'), 'w').write('\n'.join(obc) + '\n')
print(len(fes), 'FES tasks,', len(obc), 'OBC tasks')

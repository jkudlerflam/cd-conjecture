"""C160 Typhon task (Claude Opus agent, 2026-11-08): open-chain DMRG entanglement profile for the K-stratum
chiral-soliton-lattice family near h_F (PROOF section 3.4, csl_tilt.py):
    v = (-1, -1, -1 + kappa eps^2/2, 0, 0, -eps, -h eps^2, 0, 0)   (est.py coordinates: Jx Jy Jz Dx Dy Dz fx fy fz)
Usage: python3 dmrg_task.py LABEL EPS KAPPA H L OUTDIR
Writes OUTDIR/LABEL_L{L}.json with E, S(l), c_prof, c_sm (C155 dmrg_c.profile_fits), Smid, truncation history.
"""
import sys, json, time
import numpy as np
import dmrg_c as C
import dmrg as D

label, eps, kappa, h, L, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
v = [-1, -1, -1 + kappa * eps ** 2 / 2, 0, 0, -eps, -h * eps ** 2, 0, 0]
hb, edge = C.chain_bond(v)
t0 = time.time()
chis = (16, 32, 64, 96, 128, 160, 160)
r = D.dmrg(hb, L, chis=chis, seed=1, edge_field=edge)
S = r['S']
cp, cs, mid = C.profile_fits(S, L)
hist = r['hist']
rec = dict(label=label, eps=eps, kappa=kappa, h=h, L=L, v=v, E=r['E'], c_prof=cp, c_sm=cs, Smid=mid,
           S=[float(x) for x in S], hist=[[int(a), float(b), float(c)] for a, b, c in hist],
           dE_last=float(hist[-1][1] - hist[-2][1]), seconds=time.time() - t0)
json.dump(rec, open(f'{out}/{label}_L{L}.json', 'w'))
print(label, L, cp, cs, mid, hist[-1])

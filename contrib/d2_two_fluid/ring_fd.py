#!/usr/bin/env python3
"""Periodic-ring spectra and entanglement for a model spec of models.py.

usage: ring_fd.py SPEC L [--nev 6] [--out data/ring_NAME_L.json]
Records the lowest nev levels in every momentum sector k = 0..L-1 (K = 2 pi k / L), with the
mean number of up spins <N> of each level (d = 2 only), and the von Neumann and second Renyi
entropies S(l), l = 1..L-1, of the lowest state overall.
"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ring_ed import sector_levels, ground_state_full, entropies
from models import get

ap = argparse.ArgumentParser()
ap.add_argument('spec'); ap.add_argument('L', type=int)
ap.add_argument('--nev', type=int, default=6)
ap.add_argument('--out', default=None)
a = ap.parse_args()
h, d = get(a.spec)
L = a.L
t0 = time.time()
lev, Ns = {}, {}
for k in range(L):
    w, v, reps, per = sector_levels(h, L, k, a.nev, d, return_vecs=True)
    lev[k] = w.tolist()
    if d == 2:
        pop = np.array([bin(int(r)).count('1') for r in reps])
        Ns[k] = ((np.abs(v) ** 2) * pop[:, None]).sum(0).tolist()
    print(k, np.round(w, 6).tolist(), f"{time.time()-t0:.0f}s", flush=True)
k0 = min(lev, key=lambda k: lev[k][0])
E0, psi = ground_state_full(h, L, k0, d)
S1 = entropies(psi, L, d, 1); S2 = entropies(psi, L, d, 2)
out = a.out or f"data/ring_{a.spec.replace(':','_').replace(',','_')}_L{L}.json"
json.dump(dict(spec=a.spec, L=L, d=d, levels=lev, N=Ns, k0=k0, E0=E0, S1=S1.tolist(), S2=S2.tolist(),
               time_s=time.time() - t0), open(out, 'w'))
print('done', out, E0, k0, time.time() - t0)

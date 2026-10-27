#!/usr/bin/env python3
"""Ising line of the XXZ chain in a transverse field near its junction with the
Delta = -1 critical line (Dmitriev-Krivnov-Ovchinnikov, cond-mat/0112145).

H = sum_j (S^xS^x + S^yS^y + Delta S^zS^z) - h sum_j S^x  (spin-1/2), coded in Pauli
units as J = diag(1,1,Delta), field = (-2h, 0, 0) (overall factor 4).
For each Delta: locate h_c(Delta) from the crossing of L*[E0(k=pi)-E0(k=0)] for two sizes,
then measure at h_c the PBC entanglement slope (l >= L/4) and the Casimir c with the
spin-2 velocity, at L = 12..20 (even).
"""
import os, sys, json
for _k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_k, "1")
import numpy as np
from scipy.optimize import brentq
import ed

def h_of(Delta, h):
    return ed.bond(np.diag([1., 1., Delta]), [-2*h, 0, 0])

def gapL(Delta, h, L):
    hh = h_of(Delta, h)
    e0 = ed.lowest(ed.sector_hamiltonian(L, hh, 0), k=1)[0][0]
    ep = ed.lowest(ed.sector_hamiltonian(L, hh, L//2), k=1)[0][0]
    return L*(ep - e0)

def hc(Delta, L1=12, L2=16, lo=0.05, hi=3.0):
    f = lambda h: gapL(Delta, h, L1) - gapL(Delta, h, L2)
    hs = np.linspace(lo, hi, 30); fs = [f(h) for h in hs]
    for i in range(len(hs)-2, -1, -1):  # largest-h crossing: the order-disorder one
        if fs[i]*fs[i+1] < 0:
            return brentq(f, hs[i], hs[i+1], xtol=1e-6)
    return np.nan

def measure(Delta, h, Ls=(12, 14, 16, 18, 20)):
    hh = h_of(Delta, h)
    rows = []
    for L in Ls:
        secs = sorted(set([0, 1, 2, L-1, L-2, L//2, L//2+1, L//2-1, L//2+2, L//2-2]))
        sp, vec = {}, {}
        for m in secs:
            e, u = ed.lowest(ed.sector_hamiltonian(L, hh, m), k=3)
            sp[m] = e; vec[m] = u[:, 0]
        m0 = min(sp, key=lambda m: sp[m][0]); E0 = sp[m0][0]
        D2 = min(sp[m][0] for m in secs if min((m-m0) % L, (m0-m) % L) == 2) - E0
        v2 = L*D2/(4*np.pi)
        S = ed.entropies(ed.expand(L, m0, vec[m0]), L)
        ls = np.arange(int(np.ceil(L/4)), L//2+1)
        x = np.log(L/np.pi*np.sin(np.pi*ls/L))/3
        cS = np.polyfit(x, S[ls-1], 1)[0]
        rows.append(dict(L=L, E0=E0, m0=m0, v2=v2, cS=cS,
                         gap=L*(sorted(np.concatenate([sp[m] for m in secs]))[1]-E0)))
    L = np.array([r['L'] for r in rows], float); E = np.array([r['E0'] for r in rows])
    A = np.vstack([np.ones_like(L), L**-2, L**-4]).T
    e, a, b = np.linalg.lstsq(A, E/L, rcond=None)[0]
    cE = -6*a/np.pi/rows[-1]['v2']
    return rows, cE

if __name__ == '__main__':
    out = {}
    for Delta in [float(x) for x in sys.argv[1:]]:
        h = hc(Delta)
        rows, cE = measure(Delta, h)
        out[str(Delta)] = dict(h_c=h, rows=rows, c_E=cE)
        print(f"Delta={Delta:+.3f} h_c={h:.4f}  c_E={cE:.3f}  cS(L)={[round(r['cS'],3) for r in rows]}  v2={rows[-1]['v2']:.3f}  gapL={[round(r['gap'],3) for r in rows]}", flush=True)
    fn = 'results/junction.json'
    old = json.load(open(fn)) if os.path.exists(fn) else {}
    old.update(out); json.dump(old, open(fn, 'w'), indent=1)

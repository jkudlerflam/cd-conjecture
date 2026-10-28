#!/usr/bin/env python3
"""Level spectroscopy for Ising lines on rings.

For a chain with an order parameter at momentum k0 = 0 or pi, odd under the
parity P (or, for P-even order at k0 = pi, odd under one-site translation), the
"sigma" level is the lowest state of the sector (m_sig, P_sig) and the "eps"
level the first excitation of the ground-state sector (0, +1).  At an Ising
point L*(E_sig - E_0) and L*(E_eps - E_0) are size independent with ratio 1/8.
The crossing of L*(E_sig - E_0) for sizes L1 < L2 locates the line
(phenomenological renormalization).

usage: track.py FAMILY --fixed k=v ... --scan NAME LO HI N --solve NAME LO HI
                --sig 0,-1 [--L 12 16 20] [--out data/track_xxx.json]
For each value of the scan parameter, finds every root in [LO, HI] of
L1*gap(L1) - L2*gap(L2) for consecutive size pairs, and records gaps, the
ratio sigma/eps, the ring entanglement slope and the Casimir combination
at the root.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
import sys, json, argparse
import multiprocessing as mp
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ring
from scan import bond_of


def sector_levels(h, L, m, P, k=3):
    P = P if ring.has_parity(h) else None
    return ring.lowest(ring.sector_hamiltonian(L, h, m % L, P), k=k)[0]


def gaps(fam, kw, L, sig):
    """Gaps measured from the vacuum-sector ground state E0 = lowest level of (0,+):
    sigma = |E(sig) - E0| (absolute value, so that a sigma level that dips below the
    vacuum on the ordered side does not create spurious roots), eps = E1(0,+) - E0."""
    h = bond_of(fam, **kw)
    par = ring.has_parity(h)
    ms, Ps = sig
    msig = 0 if ms == 0 else L // 2
    Psig = Ps if par else 1
    vac = sector_levels(h, L, 0, 1)
    E0 = vac[0]
    s = vac if (msig, Psig) == (0, 1) else sector_levels(h, L, msig, Psig)
    esig = abs(s[1] - E0) if (msig, Psig) == (0, 1) else abs(s[0] - E0)
    eps = vac[1] - E0
    return L * esig, L * eps


def roots(fam, fixed, sname, sval, solve, L1, L2, sig, ngrid=24):
    lo, hi, yname = solve
    kw = dict(fixed); kw[sname] = sval

    def f(y):
        kw2 = dict(kw); kw2[yname] = y
        return gaps(fam, kw2, L1, sig)[0] - gaps(fam, kw2, L2, sig)[0]
    ys = np.linspace(lo, hi, ngrid)
    fs = [f(y) for y in ys]
    out = []
    for i in range(ngrid - 1):
        if np.isfinite(fs[i]) and np.isfinite(fs[i + 1]) and fs[i] * fs[i + 1] < 0:
            try:
                out.append(brentq(f, ys[i], ys[i + 1], xtol=1e-5))
            except ValueError:
                pass
    return out


def measure(fam, kw, sig, Ls):
    rows = []
    for L in Ls:
        gs, ge = gaps(fam, kw, L, sig)
        rows.append(dict(L=L, Lgap_sig=gs, Lgap_eps=ge, ratio=gs / ge))
    return rows


def work(task):
    fam, fixed, sname, sval, solve, Ls, sig = task
    res = []
    for L1, L2 in zip(Ls[:-1], Ls[1:]):
        for y in roots(fam, fixed, sname, sval, solve, L1, L2, sig):
            kw = dict(fixed); kw[sname] = sval; kw[solve[2]] = y
            res.append(dict(L1=L1, L2=L2, root=y, meas=measure(fam, kw, sig, (L1, L2))))
    return dict(scan=sval, roots=res)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fam')
    ap.add_argument('--fixed', nargs='*', default=[])
    ap.add_argument('--scan', nargs=4, required=True)
    ap.add_argument('--solve', nargs=3, required=True)
    ap.add_argument('--sig', default='0,-1')
    ap.add_argument('--L', type=int, nargs='+', default=[12, 16])
    ap.add_argument('--procs', type=int, default=os.cpu_count())
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    fixed = {s.split('=')[0]: float(s.split('=')[1]) for s in a.fixed}
    sig = tuple(int(x) for x in a.sig.split(','))
    svals = np.linspace(float(a.scan[1]), float(a.scan[2]), int(a.scan[3]))
    solve = (float(a.solve[1]), float(a.solve[2]), a.solve[0])
    for L in a.L:
        ring.orbits(L)
    tasks = [(a.fam, fixed, a.scan[0], float(s), solve, a.L, sig) for s in svals]
    with mp.get_context('fork').Pool(a.procs) as pool:
        res = pool.map(work, tasks, chunksize=1)
    json.dump(dict(fam=a.fam, fixed=fixed, scan=a.scan[0], solve=a.solve[0], sig=sig, Ls=a.L, res=res),
              open(a.out, 'w'), indent=1)
    for r in res:
        print(r['scan'], [(x['L1'], x['L2'], round(x['root'], 4), round(x['meas'][-1]['ratio'], 3)) for x in r['roots']])


if __name__ == '__main__':
    main()

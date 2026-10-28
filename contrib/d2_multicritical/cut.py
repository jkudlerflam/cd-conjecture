#!/usr/bin/env python3
"""One-dimensional cut through a family on rings: spectra, c_S and K_L.

usage: cut.py FAMILY --fixed k=v ... --x NAME LO HI N --L 12 16 20 [--k 4] --out data/cuts/X.json

Per point and size: lowest k levels of every (m, P) sector, the ground-state
sector, the chord slope c_S, and the density Luttinger estimate
K_L = 2 pi S(2pi/L) / (2pi/L) with S(q) = (1/L) sum_ij e^{iq(i-j)} <n_i n_j>_c,
n = (1 - Z)/2 (C15's estimator; free fermions give 1, the XXZ point
Delta = cos(gamma) gives pi / (2 (pi - gamma))).  K_L is only meaningful in a
c = 1 liquid whose charge is the Z magnetization (exact or emergent U(1)).
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
import sys, json, argparse
import multiprocessing as mp
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ring
from scan import bond_of


def zz_corr(psi, L):
    p = np.abs(psi) ** 2
    s = np.arange(1 << L)
    z = np.array([1 - 2 * ((s >> j) & 1) for j in range(L)], dtype=np.int8)
    zm = z.astype(float) @ p                     # <Z_j>
    C = np.array([(z[0].astype(float) * z[r]) @ p for r in range(L)])   # <Z_0 Z_r>
    return zm, C


def K_of(zm, C, L):
    n_c = (C - zm[0] * np.roll(zm, 0)[np.arange(L)]) / 4.0     # <n_0 n_r>_c (translation invariant)
    q = 2 * np.pi / L
    Sq = float(np.real(np.sum(n_c * np.exp(1j * q * np.arange(L)))))
    return 2 * np.pi * Sq / q


def work(task):
    fam, kw, L, k = task
    h = bond_of(fam, **kw)
    sp, key, psi = ring.spectrum(L, h, k=k, want_gs=True)
    S = ring.entropies(psi, L)
    zm, C = zz_corr(psi, L)
    return dict(p=kw, L=L, gs=list(key), E0=float(min(np.nanmin(v) for v in sp.values())),
                lev={f'{m},{P}': [float(x) for x in v] for (m, P), v in sp.items()},
                S=S.tolist(), cS=ring.c_chord(S, L), K=K_of(zm, C, L), mz=float(zm.mean()),
                zz=C.tolist())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fam')
    ap.add_argument('--fixed', nargs='*', default=[])
    ap.add_argument('--x', nargs=4, required=True)
    ap.add_argument('--L', type=int, nargs='+', default=[12, 16, 20])
    ap.add_argument('--k', type=int, default=4)
    ap.add_argument('--procs', type=int, default=os.cpu_count())
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    fixed = {s.split('=')[0]: float(s.split('=')[1]) for s in a.fixed}
    xs = np.linspace(float(a.x[1]), float(a.x[2]), int(a.x[3]))
    tasks = []
    for L in a.L:
        ring.orbits(L)
        for x in xs:
            kw = dict(fixed); kw[a.x[0]] = round(float(x), 10)
            tasks.append((a.fam, kw, L, a.k))
    with mp.get_context('fork').Pool(a.procs) as pool:
        res = pool.map(work, tasks, chunksize=1)
    json.dump(dict(fam=a.fam, fixed=fixed, x=a.x[0], res=res), open(a.out, 'w'))
    for r in res:
        print(r['L'], r['p'][a.x[0]], r['gs'], 'cS=%.3f K=%.3f mz=%.3f' % (r['cS'], r['K'], r['mz']), flush=True)


if __name__ == '__main__':
    main()

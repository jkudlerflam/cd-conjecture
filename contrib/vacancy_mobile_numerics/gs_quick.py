#!/usr/bin/env python3
"""Lowest levels of the best first-order candidates of one ring without momentum projection (faster than
scan.py for L >= 13).  Lanczos on the full (M, word, 2Q) sector, momenta read from <psi|T|psi>; for a
real eigenvector of a +-k pair this is cos k.  Degeneracies across momenta appear as equal eigenvalues.
usage: gs_quick.py --L 13 --J 0.25 --mu_v 0.6 [--ncand 2 --nev 6] --out f.json
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse
import numpy as np
from scipy.sparse.linalg import eigsh, LinearOperator
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'routeA_demonstrator'))
import sqed  # noqa: E402
import scan  # noqa: E402
import selector  # noqa: E402
import vacbond as VB  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--t', type=float, default=1.0)
    ap.add_argument('--J', type=float, default=0.25)
    ap.add_argument('--mu_v', type=float, default=0.6)
    ap.add_argument('--lam', type=float, default=0.15)
    ap.add_argument('--ncand', type=int, default=2)
    ap.add_argument('--nev', type=int, default=6)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    L = a.L
    ek = selector.unit_energies(a.lam, VB.DEMO['D0'], VB.DEMO['D1'], VB.DEMO['kap'], VB.DEMO['Jp'], VB.DEMO['Jz'], kmax=4)
    cands = []
    for M in range(2, L):
        ef = a.t * sqed.lowest_free(L, M, 0.0)
        for w in scan.necklaces(M):
            cands.append(dict(M=M, word=w, E1=float(ef + a.J * scan.alpha(L, M) * scan.run_energy(w, ek, VB.DEMO['mu'], a.lam) - a.mu_v * M)))
    cands.sort(key=lambda c: c['E1'])
    res = dict(args=vars(a), candidates=cands[:10], sectors=[])
    t0 = time.time()
    for c in cands[:a.ncand]:
        Q2 = c['M'] % 2
        ring = sqed.Ring(L, c['word'], Q2, a.t, a.J, a.mu_v, a.lam)
        op = LinearOperator((ring.N, ring.N), matvec=ring.H, dtype=float)
        w, U = eigsh(op, k=a.nev, which='SA', tol=1e-11, ncv=max(3 * a.nev, 30), maxiter=20000)
        o = np.argsort(w); w, U = w[o], U[:, o]
        cosk = [float(U[:, i] @ ring.T(U[:, i])) for i in range(a.nev)]
        res['sectors'].append(dict(M=c['M'], word=c['word'], Q2=Q2, dim=[ring.D, ring.C], levels=w.tolist(), cosk=cosk))
        print(f"L={L} M={c['M']} {c['word']} 2Q={Q2} dim={ring.D}x{ring.C} levels={np.round(w, 9).tolist()} "
              f"cos k={np.round(cosk, 4).tolist()} t={time.time() - t0:.0f}s", flush=True)
        json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()

"""(5) Uniform zero-port selector bound of C80 PROOF eq. (16), tested on the
complete single-cell Hilbert space (every data letter at every site, every
terminal marker; all illegal and multi-token sectors), and the marker-free
ring density bound. Uses my12.py only."""
import json, itertools, os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
from scipy.sparse.linalg import eigsh
import my12 as M

out = []
for K in (6., 10., 100.):
    for k in range(0, 6):
        p = M.default_params(max(k, 1), D=4., Delta=2.2, g=.07, td=.003, phiA=.31, musel=0., K=K)
        seeds = [w + (m,) for w in itertools.product(M.DATA, repeat=k) for m in M.MARKERS]
        bonds = [(-1, 0)] + [(i, i + 1) for i in range(k)]
        st, H, _ = M.closure(seeds, p, bonds, left='E0', cap=10 ** 6)
        if len(st) < 3000:
            E = float(np.linalg.eigvalsh(H.toarray())[0])
        else:
            E = float(eigsh(H, k=1, which='SA', tol=1e-12, return_eigenvectors=False)[0])
        Ev = -2 * np.cos(np.pi / (k + 2))
        lo = Ev - max(k - 1, 0) / (K - 4)
        r = dict(K=K, k=k, dim=len(st), E=E, lower=lo, upper=Ev, holds=bool(lo - 1e-10 <= E <= Ev + 1e-10))
        out.append(r); print(r, flush=True)
rings = []
for K in (10., 100.):
    for N in range(2, 6):
        p = M.default_params(2, D=4., Delta=2.2, g=.07, td=.003, musel=0., K=K)
        seeds = list(itertools.product(M.DATA, repeat=N))
        bonds = [(i, (i + 1) % N) for i in range(N)]
        st, H, _ = M.closure(seeds, p, bonds, cap=10 ** 6)
        E = float(eigsh(H, k=1, which='SA', tol=1e-12, return_eigenvectors=False)[0]) if len(st) > 3000 else float(np.linalg.eigvalsh(H.toarray())[0])
        r = dict(K=K, ring=N, dim=len(st), E=E, bound=-N / (K - 4), holds=bool(E >= -N / (K - 4) - 1e-10))
        rings.append(r); print(r, flush=True)
json.dump(dict(cells=out, marker_free_rings=rings), open('check_selector12.json', 'w'), indent=1)

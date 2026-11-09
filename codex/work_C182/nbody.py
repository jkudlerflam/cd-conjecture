"""C182: few-body bound states of the Fibonacci-cusp model (C173 eq. 3.2) at general (|w|, phi), in gap coordinates.
Adapted from C175's nbody_box.py (same conventions, checked against it below), with the hop table of pairband.py.
N particles x_1 < ... < x_N, gaps g_i >= 2, psi = e^{iK x_1} f(g). Box: every gap <= R (Dirichlet restriction,
so eigenvalues are upper bounds and E < threshold certifies binding).
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
from pairband import hops

def build(N, K, w, R, t2=-1.0, V2=-2.0, V3=-2/3):
    G = N-1
    shape = (R-1,)*G
    n = int(np.prod(shape))
    grid = np.indices(shape).reshape(G, n) + 2
    diag = np.zeros(n)
    for gi in range(G):
        diag += V2*(grid[gi] == 2) + V3*(grid[gi] == 3)
    rows = [np.arange(n)]; cols = [np.arange(n)]; vals = [diag.astype(complex)]
    for d, ad in hops(w, t2).items():
        for i in range(N):
            new = grid.copy(); ph = 1.0+0j
            if i == 0:
                new[0] += d; ph = np.exp(-1j*K*d)
            elif i == N-1:
                new[G-1] -= d
            else:
                new[i-1] -= d; new[i] += d
            ok = np.all((new >= 2) & (new <= R), axis=0)
            src = np.nonzero(ok)[0]
            tgt = np.ravel_multi_index(tuple(new[:, src]-2), shape)
            rows.append(src); cols.append(tgt); vals.append(np.full(len(src), ad*ph))
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    assert abs(H-H.getH()).max() < 1e-10
    return H

def lowest(N, K, w, R, k=1, v0=None, **kw):
    H = build(N, K, w, R, **kw)
    e, v = spla.eigsh(H, k=k, which='SA', tol=1e-11, v0=v0)
    o = np.argsort(e)
    return e[o], v[:, o]

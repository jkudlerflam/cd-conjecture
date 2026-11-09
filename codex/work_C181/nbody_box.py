"""C175: few-body bound states of model (3.2) on the infinite line, in relative (gap) coordinates.
N particles x_1 < ... < x_N, gaps g_i = x_{i+1} - x_i >= 2 (blockade), psi = e^{iK x_1} f(g_1..g_{N-1}).
Particle i moved by -d (amplitude a_d, i.e. psi(.., x_i - d, ..)):
  i = 1:      g_1 -> g_1 + d, phase e^{-iKd}
  1 < i < N:  g_{i-1} -> g_{i-1} - d, g_i -> g_i + d
  i = N:      g_{N-1} -> g_{N-1} - d
Diagonal: V2 for each gap = 2, V3 for each gap = 3 (pairs at distance 2 or 3 are always neighbours in order,
since gaps >= 2).  Box: each gap <= R (hard wall far away), so bound states are exact up to e^{-R/size}.
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, itertools, sys
from scat import valleys, hops

def build(N, K, w, R, t2=-1.0, V2=-2.0, V3=-2/3):
    G = N-1
    shape = (R-1,)*G               # gap values 2..R
    n = int(np.prod(shape))
    grid = np.indices(shape).reshape(G, n) + 2   # gaps
    diag = np.zeros(n)
    for gi in range(G):
        diag += V2*(grid[gi] == 2) + V3*(grid[gi] == 3)
    rows = [np.arange(n)]; cols = [np.arange(n)]; vals = [diag.astype(complex)]
    a = hops(w, t2)
    for d, ad in a.items():
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
            # (H f)(src) += a_d * ph * f(tgt)
            rows.append(src); cols.append(tgt); vals.append(np.full(len(src), ad*ph))
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    herm = abs(H-H.getH()).max()
    assert herm < 1e-10, herm
    return H

def lowest(N, K, w, R, k=1, **kw):
    H = build(N, K, w, R, **kw)
    e = spla.eigsh(H, k=k, which='SA', tol=1e-10)[0]
    return np.sort(e)

if __name__ == "__main__":
    wa = float(sys.argv[1]); N = int(sys.argv[2]); R = int(sys.argv[3])
    Ks = [float(x) for x in sys.argv[4].split(',')]
    w = 1j*wa
    (ka, ea), (kb, eb) = valleys(w)[:2]
    for K in Ks:
        e = lowest(N, K, w, R)
        print(f"|w|={wa} N={N} R={R} K={K:.4f}: E0 = {e[0]:.6f}, E0 - N eps_min = {e[0]-N*ea:.6f}", flush=True)

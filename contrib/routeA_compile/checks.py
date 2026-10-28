"""Exact-diagonalization checks of the dimer-port cell (see README).

For each (n, cells, K) it verifies
  * the code (all dimers in 'none') is an exact zero-energy kernel of H0
    with a gap, and G^dagger V G = 0;
  * the second-order bond H2 = -P V R V P (exact resolvent, truncated
    space exact at this order) and the cubic coefficient H3;
  * the des Cloizeaux effective Hamiltonian from the lowest eigenvectors
    of H0 + lambda V, and the residual against lambda^2 H2 (should scale
    as lambda^3 and match lambda^3 H3);
  * Gaussianity of H2: its one-body part (from the 0- and 1-fermion
    sectors) reproduces it up to a non-Gaussian remainder.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import json
import sys
import time
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from model import Chain


def second_third(ch, H0, V):
    code = ch.code
    mask = np.ones(ch.dim, bool)
    mask[code] = False
    Q = np.where(mask)[0]
    VQP = V[Q][:, code].toarray()
    HQQ = H0[Q][:, Q].tocsc()
    X = sla.spsolve(HQQ, VQP)
    if X.ndim == 1:
        X = X[:, None]
    H2 = -VQP.conj().T @ X
    VQQ = V[Q][:, Q]
    H3 = X.conj().T @ (VQQ @ X)
    return 0.5 * (H2 + H2.conj().T), 0.5 * (H3 + H3.conj().T)


def gap(ch, H0, nev=6):
    """Lowest H0 eigenvalue outside the code (carrier sectors)."""
    mask = np.ones(ch.dim, bool)
    mask[ch.code] = False
    Q = np.where(mask)[0]
    HQQ = H0[Q][:, Q]
    if HQQ.shape[0] < 4000:
        return float(la.eigvalsh(HQQ.toarray())[0])
    return float(sla.eigsh(HQQ, k=1, which='SA', return_eigenvectors=False)[0])


def effective(ch, H0, V, lam):
    nc = len(ch.code)
    Heff = np.zeros((nc, nc), complex)
    pos = {int(c): q for q, c in enumerate(ch.code)}
    H = (H0 + lam * V).tocsr()
    nxt = np.inf
    for N in np.unique(ch.number[ch.code]):
        idx = np.where(ch.number == N)[0]
        cidx = [c for c in ch.code if ch.number[c] == N]
        k = len(cidx)
        Hs = H[idx][:, idx]
        if len(idx) <= 3000:
            E, U = la.eigh(Hs.toarray())
        else:
            E, U = sla.eigsh(Hs, k=k + 2, which='SA', tol=1e-13)
            o = np.argsort(E)
            E, U = E[o], U[:, o]
        nxt = min(nxt, E[k])
        E, U = E[:k], U[:, :k]
        loc = {int(g): r for r, g in enumerate(idx)}
        rows = [loc[int(c)] for c in cidx]
        O = U[rows, :]
        s, w = la.eigh(O.conj().T @ O)
        W = O @ (w @ np.diag(s ** -0.5) @ w.conj().T)
        h = W @ np.diag(E) @ W.conj().T
        q = [pos[int(c)] for c in cidx]
        Heff[np.ix_(q, q)] = h
    return Heff, nxt


def one_body(ch, H2):
    """E0 and the one-body matrix h_ab of H2 in the physical data operators."""
    dm = ch.data_modes
    vac = ch.code_op([])  # identity
    nc = len(ch.code)
    # vacuum code state = no data fermions
    vacidx = [q for q, c in enumerate(ch.code)
              if all(((ch.states[c] >> m) & 1) == 0 for m in dm)][0]
    e = np.zeros(nc)
    e[vacidx] = 1
    E0 = float(np.real(e @ H2 @ e))
    kets = [ch.code_op([(m, 1)]) @ e for m in dm]
    h = np.array([[k1.conj() @ H2 @ k2 for k2 in kets] for k1 in kets])
    h -= E0 * np.eye(len(dm))
    G = E0 * np.eye(nc, dtype=complex)
    for a, ma in enumerate(dm):
        for b, mb in enumerate(dm):
            if abs(h[a, b]) > 0:
                G += h[a, b] * ch.code_op([(ma, 1), (mb, 0)])
    return E0, h, G


def run(n, M, K, lams, params=None, verbose=True):
    t0 = time.time()
    ch = Chain(n, M, K, params)
    H0, V = ch.hamiltonians()
    out = dict(n=n, cells=M, K=K, dim=ch.dim, code_dim=len(ch.code),
               params=ch.p)
    out['herm_H0'] = float(abs(H0 - H0.getH()).max())
    out['herm_V'] = float(abs(V - V.getH()).max())
    out['code_residual'] = float(abs(H0[:, ch.code]).max()) if H0[:, ch.code].nnz else 0.0
    out['first_order'] = float(abs(V[ch.code][:, ch.code]).max()) if V[ch.code][:, ch.code].nnz else 0.0
    out['gap'] = gap(ch, H0)
    H2, H3 = second_third(ch, H0, V)
    E0, h, G = one_body(ch, H2)
    tr = np.trace(H2).real / len(ch.code)
    out['H2_norm_traceless'] = float(la.norm(H2 - tr * np.eye(len(ch.code)), 2))
    out['nonGaussian_H2'] = float(la.norm(H2 - G, 2))
    out['H3_norm'] = float(la.norm(H3, 2))
    out['one_body_h'] = [[complex(z).real, complex(z).imag] for z in h.ravel()]
    rows = []
    for lam in lams:
        Heff, nxt = effective(ch, H0, V, lam)
        E = la.norm(Heff - lam ** 2 * H2, 2)
        E3 = la.norm(Heff - lam ** 2 * H2 - lam ** 3 * H3, 2)
        d = Heff - lam ** 2 * H2
        Etr = la.norm(d - np.trace(d).real / len(d) * np.eye(len(d)), 2)
        rows.append(dict(lam=lam, E=float(E), E_over_lam3=float(E / lam ** 3),
                         E_traceless=float(Etr), E_after_H3=float(E3),
                         E_after_H3_over_lam4=float(E3 / lam ** 4),
                         next_level=float(nxt)))
    out['ed'] = rows
    out['seconds'] = time.time() - t0
    if verbose:
        print(json.dumps({k: v for k, v in out.items() if k not in ('one_body_h', 'params')},
                         indent=1), flush=True)
    return out, ch, H2, h


if __name__ == '__main__':
    n, M, K = map(int, sys.argv[1:4])
    run(n, M, K, [0.08, 0.04, 0.02])

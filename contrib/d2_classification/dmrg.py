#!/usr/bin/env python3
"""Minimal two-site finite DMRG (open chain) for a general nearest-neighbour qubit chain.

H = sum_{j<L-1} sum_ab J_ab s^a_j s^b_{j+1} + sum_j h.s_j  (Pauli, open boundaries,
uniform on-site field).  No symmetries are used, so DM terms, oblique fields and
incommensurate spirals are all allowed.  Output: ground energy, truncation error, and the
entanglement entropy S(l) at every bond, from the last left-to-right sweep.

usage: dmrg.py NAME --L 64 --chi 32 64 128 192 --sweeps 2 [--seed 1] [--out results/dmrg_NAME_L64.json]
"""
import os
for _k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_k, "4")
import argparse, json, time
import numpy as np
from scipy.sparse.linalg import eigsh, LinearOperator

PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)


def mpo(J, h):
    J = np.asarray(J, float); h = np.asarray(h, float)
    D = 5
    W = np.zeros((D, D, 2, 2), complex)
    I = np.eye(2)
    W[0, 0] = I
    W[4, 4] = I
    W[4, 0] = np.einsum('a,aij->ij', h, PAULI)
    for a in range(3):
        W[4, a + 1] = PAULI[a]
        W[a + 1, 0] = np.einsum('b,bij->ij', J[a], PAULI)
    return W


def bond_op(J, h):
    J = np.asarray(J, float); h = np.asarray(h, float)
    I = np.eye(2)
    B = np.zeros((4, 4), complex)
    for a in range(3):
        B += h[a] * 0.5 * (np.kron(PAULI[a], I) + np.kron(I, PAULI[a]))
        for b in range(3):
            B += J[a, b] * np.kron(PAULI[a], PAULI[b])
    return B


def apply2(Lenv, W, Renv, th):
    t = np.tensordot(Lenv, th, axes=([2], [0]))          # a' w s1 s2 b
    t = np.tensordot(t, W, axes=([1, 2], [0, 3]))         # a' s2 b u s1'
    t = np.tensordot(t, W, axes=([3, 1], [0, 3]))         # a' b s1' v s2'
    t = np.tensordot(t, Renv, axes=([1, 3], [2, 1]))      # a' s1' s2' b'
    return t


def grow_left(Lenv, W, A):
    t = np.tensordot(Lenv, A, axes=([2], [0]))            # a' w s b
    t = np.tensordot(t, W, axes=([1, 2], [0, 3]))         # a' b w' s'
    t = np.tensordot(A.conj(), t, axes=([0, 1], [0, 3]))  # b' b w'
    return t.transpose(0, 2, 1)


def grow_right(Renv, W, B):
    t = np.tensordot(B, Renv, axes=([2], [2]))            # a s b' w'
    t = np.tensordot(t, W, axes=([1, 3], [3, 1]))         # a b' w s'
    t = np.tensordot(B.conj(), t, axes=([1, 2], [3, 1]))  # a' a w
    return t.transpose(0, 2, 1)


def dmrg(J, h, L, chis, sweeps=2, seed=1, tol=1e-9, cutoff=1e-12, verbose=True):
    rng = np.random.default_rng(seed)
    W = mpo(J, h)
    D = W.shape[0]
    # random right-canonical MPS with bond 8
    chi0 = 8
    dims = [1]
    for i in range(1, L):
        dims.append(min(chi0, 2 ** i, 2 ** (L - i)))
    dims.append(1)
    M = [rng.normal(size=(dims[i], 2, dims[i + 1])) + 1j * rng.normal(size=(dims[i], 2, dims[i + 1])) for i in range(L)]
    for i in range(L - 1, 0, -1):
        a, s, b = M[i].shape
        q, r = np.linalg.qr(M[i].reshape(a, s * b).T)
        M[i] = q.T.reshape(-1, s, b)
        M[i - 1] = np.tensordot(M[i - 1], r.T, axes=([2], [0]))
    M[0] /= np.linalg.norm(M[0])
    Lenv = [None] * (L + 1)
    Renv = [None] * (L + 1)
    Lenv[0] = np.zeros((1, D, 1), complex); Lenv[0][0, 4, 0] = 1
    Renv[L] = np.zeros((1, D, 1), complex); Renv[L][0, 0, 0] = 1
    for i in range(L - 1, 0, -1):
        Renv[i] = grow_right(Renv[i + 1], W, M[i])
    E = None
    hist = []
    S = np.zeros(L - 1)
    Eb = np.zeros(L - 1)
    B2 = bond_op(J, h).reshape(2, 2, 2, 2)
    for chi in chis:
        for sw in range(sweeps):
            maxtr = 0.0
            # left -> right
            for i in range(L - 1):
                th = np.tensordot(M[i], M[i + 1], axes=([2], [0]))
                shp = th.shape
                n = th.size
                if n <= 64:
                    Hd = np.zeros((n, n), complex)
                    for k in range(n):
                        v = np.zeros(n, complex); v[k] = 1
                        Hd[:, k] = apply2(Lenv[i], W, Renv[i + 2], v.reshape(shp)).ravel()
                    w, U = np.linalg.eigh((Hd + Hd.conj().T) / 2)
                    E, v0 = w[0], U[:, 0]
                else:
                    op = LinearOperator((n, n), matvec=lambda x: apply2(Lenv[i], W, Renv[i + 2], x.reshape(shp)).ravel(), dtype=complex)
                    w, U = eigsh(op, k=1, which='SA', v0=th.ravel(), tol=tol, ncv=min(n - 1, 24))
                    E, v0 = w[0], U[:, 0]
                th = v0.reshape(shp[0] * 2, 2 * shp[3])
                u, s, vh = np.linalg.svd(th, full_matrices=False)
                keep = min(chi, max(1, int((s > cutoff * s[0]).sum())))
                tr = float((s[keep:] ** 2).sum())
                maxtr = max(maxtr, tr)
                s = s[:keep] / np.linalg.norm(s[:keep])
                p = s ** 2
                S[i] = float(-(p[p > 1e-16] * np.log(p[p > 1e-16])).sum())
                M[i] = u[:, :keep].reshape(shp[0], 2, keep)
                M[i + 1] = (s[:, None] * vh[:keep]).reshape(keep, 2, shp[3])
                Lenv[i + 1] = grow_left(Lenv[i], W, M[i])
            Sright = S.copy()
            # right -> left
            for i in range(L - 2, -1, -1):
                th = np.tensordot(M[i], M[i + 1], axes=([2], [0]))
                shp = th.shape
                n = th.size
                if n <= 64:
                    Hd = np.zeros((n, n), complex)
                    for k in range(n):
                        v = np.zeros(n, complex); v[k] = 1
                        Hd[:, k] = apply2(Lenv[i], W, Renv[i + 2], v.reshape(shp)).ravel()
                    w, U = np.linalg.eigh((Hd + Hd.conj().T) / 2)
                    E, v0 = w[0], U[:, 0]
                else:
                    op = LinearOperator((n, n), matvec=lambda x: apply2(Lenv[i], W, Renv[i + 2], x.reshape(shp)).ravel(), dtype=complex)
                    w, U = eigsh(op, k=1, which='SA', v0=th.ravel(), tol=tol, ncv=min(n - 1, 24))
                    E, v0 = w[0], U[:, 0]
                th4 = v0.reshape(shp)
                Eb[i] = float(np.real(np.vdot(th4, np.einsum('stuv,auvb->astb', B2, th4))))
                th = v0.reshape(shp[0] * 2, 2 * shp[3])
                u, s, vh = np.linalg.svd(th, full_matrices=False)
                keep = min(chi, max(1, int((s > cutoff * s[0]).sum())))
                tr = float((s[keep:] ** 2).sum())
                maxtr = max(maxtr, tr)
                s = s[:keep] / np.linalg.norm(s[:keep])
                p = s ** 2
                S[i] = float(-(p[p > 1e-16] * np.log(p[p > 1e-16])).sum())
                M[i + 1] = vh[:keep].reshape(keep, 2, shp[3])
                M[i] = (u[:, :keep] * s[None, :]).reshape(shp[0], 2, keep)
                Renv[i + 1] = grow_right(Renv[i + 2], W, M[i + 1])
            hist.append(dict(chi=chi, sweep=sw, E=float(E), maxtrunc=maxtr))
            if verbose:
                print(f"chi={chi} sweep={sw} E={E:.12f} E/L={E/L:.10f} trunc={maxtr:.2e} Smid={S[L//2-1]:.5f}", flush=True)
    return dict(E=float(E), S=S.tolist(), Eb=Eb.tolist(), hist=hist)


def fit_obc(S, L, frac=0.125, Eb=None):
    """S(l) = (c/6) ln[(2L/pi) sin(pi l/L)] + b [+ beta (Eb(l) - mean)] over l in [frac L, (1-frac) L].
    The bond-energy term (Laflorencie et al. 2006) absorbs boundary-induced oscillations."""
    S = np.asarray(S)
    l = np.arange(1, L)
    sel = (l >= frac * L) & (l <= (1 - frac) * L)
    x = np.log(2 * L / np.pi * np.sin(np.pi * l[sel] / L)) / 6
    cols = [x, np.ones_like(x)]
    if Eb is not None:
        e = np.asarray(Eb)[sel]
        cols.append(e - e.mean())
    A = np.vstack(cols).T
    return float(np.linalg.lstsq(A, S[sel], rcond=None)[0][0])


if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from controls import local_term, CONTROLS
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--chi', type=int, nargs='+', default=[32, 64, 128])
    ap.add_argument('--sweeps', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    import json as _j
    pars = _j.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'candidates_params.json')))
    if a.name in pars:
        J, h = pars[a.name]['J'], pars[a.name]['field']
    elif a.name in CONTROLS:
        J, h = CONTROLS[a.name][0], CONTROLS[a.name][1]
    else:
        ex = _j.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'extra_points.json')))
        J, h = ex[a.name]['J'], ex[a.name]['field']
    t0 = time.time()
    r = dmrg(np.asarray(J, float), np.asarray(h, float), a.L, a.chi, a.sweeps, a.seed)
    r.update(name=a.name, L=a.L, chis=a.chi, seed=a.seed, time_s=time.time() - t0,
             c_obc={str(f): fit_obc(r['S'], a.L, f) for f in (0.125, 0.25)},
             c_obc_Eb={str(f): fit_obc(r['S'], a.L, f, r['Eb']) for f in (0.125, 0.25)})
    out = a.out or f"results/dmrg_{a.name}_L{a.L}.json"
    _j.dump(r, open(out, 'w'))
    print('c_obc', r['c_obc'], 'c_obc_Eb', r['c_obc_Eb'], f"time {r['time_s']:.0f}s")

#!/usr/bin/env python3
"""Ring spectra of a logical chain by exact diagonalization (Lanczos, matrix free).
usage: ring_ed.py SPEC --Ls 4 5 6 --dN -1 0 1 [--k 6] [--out file.json]
H = sum_{j=0}^{L-1} K_{j, j+1 mod L} as a spin chain of local dimension d = 2^n (so the fermion
ring has the twist pi (N - 1), N the total logical fermion number).  For each L the sectors
N = N0 + dN are diagonalized, N0 = round(L * total filling of the Hartree-Fock chain); the
lowest k levels and the ground-state block entropies S(l), l = 1..L-1, are written.
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28."""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh
import models34 as Mo


def apply_H(psi, K4, L, d):
    """psi shape (d,)*L; K4 = K.reshape(d, d, d, d) acting on (j, j+1)."""
    out = np.zeros_like(psi)
    for j in range(L):
        j2 = (j + 1) % L
        t = np.tensordot(K4, psi, ([2, 3], [j, j2]))       # (a', b', rest...)
        # move the two new axes back to positions j, j2
        rest = [x for x in range(L) if x not in (j, j2)]
        order = [None] * L
        order[j] = 0; order[j2] = 1
        for i, x in enumerate(rest):
            order[x] = 2 + i
        out += np.transpose(t, order)
    return out


def number_mask(L, n):
    d = 2 ** n
    pc = np.array([bin(x).count('1') for x in range(d)])
    N = np.zeros((d,) * L, dtype=np.int16)
    for j in range(L):
        shape = [1] * L; shape[j] = d
        N = N + pc.reshape(shape)
    return N


def block_entropy(v, L, d, l):
    M = v.reshape(d ** l, d ** (L - l))
    s = np.linalg.svd(M, compute_uv=False)
    p = s ** 2; p = p[p > 1e-16]; p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec')
    ap.add_argument('--Ls', type=int, nargs='+', default=[4, 5, 6])
    ap.add_argument('--dN', type=int, nargs='+', default=[-1, 0, 1])
    ap.add_argument('--k', type=int, default=6)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    K, info = Mo.load(a.spec)
    n = info['n']; d = 2 ** n
    K4 = K.reshape(d, d, d, d)
    # filling of the Hartree-Fock chain (sum over sorted bands)
    import analyze34 as A
    lev = info['level']
    if lev.startswith('p'):
        _, Tp, Mp, _ = Mo.proxy(':'.join(a.spec.split(':')[:5]))
        T, M = Tp, Mp
    else:
        _, _, T, M, _ = Mo.base(':'.join(a.spec.split(':')[:5]))
    nu = sum(b['filling'] for b in A.band_data(T, M))
    tag = a.spec.replace(':', '_').replace(',', '-')
    out = a.out or os.path.join(Mo.HERE, 'data', 'ring', f'ring_{tag}.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(info=info, filling=nu, results=[])
    for L in a.Ls:
        Nmask = number_mask(L, n)
        N0 = int(round(nu * L))
        for dN in a.dN:
            N = N0 + dN
            sel = (Nmask == N)
            dim = int(sel.sum())
            if dim == 0:
                continue
            t0 = time.time()
            shape = (d,) * L

            def mv(x):
                psi = np.zeros(shape, complex)
                psi[sel] = x
                return apply_H(psi, K4, L, d)[sel]
            op = LinearOperator((dim, dim), matvec=mv, dtype=complex)
            k = min(a.k, dim - 2) if dim > 3 else 1
            if dim <= 400:
                Hm = np.array([mv(e) for e in np.eye(dim, dtype=complex)]).T
                w, U = np.linalg.eigh(0.5 * (Hm + Hm.conj().T))
                w = w[:a.k]; U = U[:, :a.k]
            else:
                v0 = np.random.default_rng(1).normal(size=dim) + 0j
                w, U = eigsh(op, k=k, which='SA', v0=v0, tol=1e-12, ncv=max(2 * k + 1, 24), maxiter=5000)
                o = np.argsort(w); w = w[o]; U = U[:, o]
            psi = np.zeros(shape, complex); psi[sel] = U[:, 0]
            S = [block_entropy(psi, L, d, l) for l in range(1, L)]
            rec = dict(L=L, N=N, dim=dim, levels=[float(x) for x in w], S=S, time_s=time.time() - t0)
            data['results'].append(rec)
            print(f"L={L} N={N} dim={dim} E0={w[0]:.12f} gaps={np.round(w[1:] - w[0], 8).tolist()} S={np.round(S, 4).tolist()} t={rec['time_s']:.0f}s", flush=True)
            json.dump(data, open(out, 'w'))


if __name__ == '__main__':
    main()

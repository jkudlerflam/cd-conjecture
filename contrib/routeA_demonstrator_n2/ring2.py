#!/usr/bin/env python3
"""Ring spectra of the two-flavour chain, resolved by both U(1) charges (N1, N2) and by momentum.

micro: ring of M cells (P_L B P_R, 324 states each) = the full 13-letter chain on L = 3M sites in its period-3
       type sector, H = sum_j H0_j + lam sum_j V(x_R of j, x_L of j+1) (the -mu per B is dropped), no truncation.
       Block-sparse storage over the cell charges, as in ../routeA_demonstrator/ring_ed.py.
eff  : ring of M logical cells (d = 4) with the exact two-cell bond at the same lam (all orders within two cells).
usage: ring2.py --M 3 --lam 0.15 --N 2 2 --m 0 1 --nev 3 [--model eff] [--out f.json]
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse, itertools
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_k, '8')
import numpy as np
from scipy.sparse.linalg import eigsh, LinearOperator
import scipy.sparse as sps
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cell2 as C

PAR = dict(DN=1.0, DC=1.0, kap=(0.6, 0.6), sg=-1.0, DCd=(0.0, -0.2), Jp=1.0, Jz=1.0)


class Cell:
    def __init__(self, H0, qc, terms):
        """terms: list of (coef, A (on cell j), B (on cell j+1), dA, dB) with charge steps as 2-vectors."""
        labels = sorted(set(map(tuple, np.round(qc).astype(int))))
        self.labels = labels
        self.lab = {l: i for i, l in enumerate(labels)}
        q = [tuple(x) for x in np.round(qc).astype(int)]
        self.idx = {l: np.array([i for i in range(len(q)) if q[i] == l]) for l in labels}
        self.dim = {l: len(self.idx[l]) for l in labels}
        self.H0 = {l: H0[np.ix_(self.idx[l], self.idx[l])] for l in labels}
        self.terms = []
        for cf, A, B, dA, dB in terms:
            Ab, Bb = {}, {}
            for l in labels:
                la = (l[0] + dA[0], l[1] + dA[1]); lb = (l[0] + dB[0], l[1] + dB[1])
                if la in self.idx:
                    m = A[np.ix_(self.idx[la], self.idx[l])]
                    if np.abs(m).max() > 1e-14: Ab[l] = (la, m)
                if lb in self.idx:
                    m = B[np.ix_(self.idx[lb], self.idx[l])]
                    if np.abs(m).max() > 1e-14: Bb[l] = (lb, m)
            self.terms.append((cf, Ab, Bb))


def micro_cell(lam):
    H0, ports, code, qc = C.cell_ops(PAR['DN'], PAR['DC'], PAR['kap'], PAR['sg'], PAR['DCd'], 0.0)
    steps = {'up1': (1, 0), 'dn1': (-1, 0), 'up2': (0, 1), 'dn2': (0, -1), 'D1': (0, 0), 'D2': (0, 0)}
    terms = [(lam * cf, ports['R'][r], ports['L'][l], steps[r], steps[l]) for cf, r, l in C.bond_terms(PAR['Jp'], PAR['Jz'])]
    return Cell(H0, qc, terms)


class Ring:
    def __init__(self, cell, M, N):
        self.c, self.M = cell, M
        self.keys = [k for k in itertools.product(cell.labels, repeat=M)
                     if (sum(x[0] for x in k), sum(x[1] for x in k)) == tuple(N)]
        self.shapes = {k: tuple(cell.dim[t] for t in k) for k in self.keys}
        self.offs, o = {}, 0
        for k in self.keys:
            self.offs[k] = o; o += int(np.prod(self.shapes[k]))
        self.N = o

    def unpack(self, v):
        return {k: v[self.offs[k]:self.offs[k] + int(np.prod(self.shapes[k]))].reshape(self.shapes[k]) for k in self.keys}

    def apply_H(self, v):
        c, M = self.c, self.M
        X = self.unpack(v)
        out = np.zeros(self.N, complex)
        Y = self.unpack(out)
        for k, x in X.items():
            y = Y[k]
            for j in range(M):
                y += np.moveaxis(np.tensordot(c.H0[k[j]], x, ([1], [j])), 0, j)
            for j in range(M):
                jn = (j + 1) % M
                for cf, Ab, Bb in c.terms:
                    if k[j] not in Ab or k[jn] not in Bb:
                        continue
                    la, A = Ab[k[j]]; lb, B = Bb[k[jn]]
                    nk = list(k); nk[j] = la; nk[jn] = lb; nk = tuple(nk)
                    if nk not in Y:
                        continue
                    t = np.moveaxis(np.tensordot(A, x, ([1], [j])), 0, j)
                    t = np.moveaxis(np.tensordot(B, t, ([1], [jn])), 0, jn)
                    Y[nk] += cf * t
        return out

    def translate(self, v):
        X = self.unpack(v)
        out = np.zeros(self.N, complex)
        Y = self.unpack(out)
        for k, x in X.items():
            nk = tuple(k[(i - 1) % self.M] for i in range(self.M))
            Y[nk][...] = np.moveaxis(x, list(range(self.M)), [(i + 1) % self.M for i in range(self.M)])
        return out

    def proj(self, v, m):
        """Exact momentum projector P_k, k = 2 pi m / M (complex)."""
        k = 2 * np.pi * m / self.M
        out = np.zeros_like(v); w = v.copy()
        for r in range(self.M):
            out += np.exp(-1j * k * r) * w
            w = self.translate(w)
        return out / self.M


def lowest(apply, proj, N, nev, shift=10.0, seed=0):
    rng = np.random.default_rng(seed)
    v0 = proj(rng.normal(size=N) + 1j * rng.normal(size=N))
    op = LinearOperator((N, N), matvec=lambda v: proj(apply(v) - shift * v), dtype=complex)
    w = eigsh(op, k=nev, which='SA', v0=v0, tol=1e-10, ncv=max(2 * nev + 10, 24), maxiter=20000, return_eigenvectors=False)
    return np.sort(w.real) + shift


def eff_ring(M, N, m, Heff, nev):
    """Dense-bond logical ring (d = 4 per cell, basis (n1, n2)), sector (N1, N2), momentum m."""
    d = 4
    st = [s for s in itertools.product(range(d), repeat=M) if (sum(x >> 1 for x in s), sum(x & 1 for x in s)) == tuple(N)]
    ix = {s: i for i, s in enumerate(st)}
    H = Heff.reshape(d, d, d, d)
    rows, cols, vals = [], [], []
    for i, s in enumerate(st):
        for j in range(M):
            jn = (j + 1) % M
            a, b = s[j], s[jn]
            for a2 in range(d):
                for b2 in range(d):
                    val = H[a2, b2, a, b]
                    if abs(val) < 1e-15: continue
                    t = list(s); t[j] = a2; t[jn] = b2; t = tuple(t)
                    if t in ix:
                        rows.append(ix[t]); cols.append(i); vals.append(val)
    Hs = sps.csr_matrix((vals, (rows, cols)), shape=(len(st), len(st)))
    T = sps.csr_matrix((np.ones(len(st)), ([ix[tuple(s[(i - 1) % M] for i in range(M))] for s in st], range(len(st)))),
                       shape=(len(st), len(st)))
    k = 2 * np.pi * m / M
    Pm = sum(np.exp(-1j * k * r) * (T ** r) for r in range(M)) / M
    A = (Pm @ Hs @ Pm).toarray() if len(st) < 3000 else None
    if A is not None:
        P = Pm.toarray()
        # eigenvalues restricted to the range of P
        U, s, _ = np.linalg.svd(P)
        Ur = U[:, s > 0.5]
        w = np.linalg.eigvalsh(Ur.conj().T @ Hs.toarray() @ Ur)
        return np.sort(w)[:nev]
    op = LinearOperator(Hs.shape, matvec=lambda v: Pm @ (Hs @ v - 10 * v), dtype=complex)
    v0 = Pm @ (np.random.default_rng(0).normal(size=len(st)) + 0j)
    w = eigsh(op, k=nev, which='SA', v0=v0, tol=1e-12, return_eigenvectors=False)
    return np.sort(w.real) + 10


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='micro')
    ap.add_argument('--M', type=int, default=3)
    ap.add_argument('--lam', type=float, default=0.15)
    ap.add_argument('--N', type=int, nargs='+', action='append', default=None, help='sector N1 N2 (repeatable)')
    ap.add_argument('--m', type=int, nargs='+', default=None)
    ap.add_argument('--nev', type=int, default=3)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    M = a.M
    sectors = a.N or [[M // 2, (M + 1) // 2]]
    ms = a.m if a.m is not None else list(range(M))
    res = dict(args=vars(a), levels={})
    t0 = time.time()
    if a.model == 'micro':
        cell = micro_cell(a.lam)
    else:
        import effective2 as E
        Heff, _, _ = E.exact_eff(dict(E.DEFAULT), a.lam)
    for N in sectors:
        if a.model == 'micro':
            ring = Ring(cell, M, N)
        for m in ms:
            if a.model == 'micro':
                w = lowest(ring.apply_H, lambda v: ring.proj(v, m), ring.N, a.nev)
                dim = ring.N
            else:
                w = eff_ring(M, N, m, Heff, a.nev); dim = None
            res['levels'][f'{N[0]}_{N[1]}_{m}'] = w.tolist()
            print(f'{a.model} M={M} N={N} m={m} dim={dim} levels={np.round(w, 10).tolist()} t={time.time() - t0:.0f}s', flush=True)
            if a.out:
                json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()

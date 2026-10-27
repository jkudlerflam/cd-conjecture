#!/usr/bin/env python3
"""Momentum-resolved exact diagonalization for translation-invariant qubit rings.

The Hamiltonian is H = sum_j T^j h T^-j with h an r-site operator on sites
(j, ..., j+r-1) mod L (r = 2 for the nearest-neighbour chains of the project,
r = 3 for the cluster-Ising control).  Site j is bit j of the basis label, and
inside the 2^r x 2^r matrix h, site j is the most significant factor.

Written for contrib/d2_classification (Claude Opus 5.5, 2026-10-27).  The orbit
construction follows the same logic as codex/work_C4_D2/d2.py but the code is
independent and is checked against full-space diagonalization in selftest().
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
from functools import lru_cache
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from numba import njit

PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)
I2 = np.eye(2)


def bond(J, field):
    """h = sum_ab J_ab s^a (x) s^b + sum_a h_a (s^a (x) 1 + 1 (x) s^a)/2 (Pauli)."""
    J = np.asarray(J, float); field = np.asarray(field, float)
    h = np.zeros((4, 4), complex)
    for a in range(3):
        h += field[a] * 0.5 * (np.kron(PAULI[a], I2) + np.kron(I2, PAULI[a]))
        for b in range(3):
            h += J[a, b] * np.kron(PAULI[a], PAULI[b])
    return h


@njit(cache=True)
def _orbits(L):
    N = 1 << L
    oid = np.full(N, -1, np.int32)
    shift = np.zeros(N, np.int8)
    reps = np.empty(N // L + 64 if L > 0 else 1, np.int64)
    reps = np.empty(N, np.int64)
    per = np.empty(N, np.int16)
    nr = 0
    mask = N - 1
    for s in range(N):
        if oid[s] >= 0:
            continue
        t = s
        p = 0
        while True:
            oid[t] = nr
            shift[t] = p
            p += 1
            t = ((t << 1) & mask) | (t >> (L - 1))
            if t == s:
                break
        reps[nr] = s
        per[nr] = p
        nr += 1
    return reps[:nr].copy(), per[:nr].copy(), oid, shift


@lru_cache(maxsize=4)
def orbits(L):
    return _orbits(L)


def sector(L, m):
    reps, per, oid, shift = orbits(L)
    allowed = (per.astype(np.int64) * m) % L == 0
    active = np.flatnonzero(allowed).astype(np.int64)
    inv = np.full(len(reps), -1, np.int64)
    inv[active] = np.arange(len(active))
    return active, inv


def _local_table(h, tol=1e-14):
    """For each input local configuration a: list of (b, h[b,a]) with h[b,a] != 0."""
    D = h.shape[0]
    cnt = np.zeros(D, np.int64)
    outs = np.zeros((D, D), np.int64)
    vals = np.zeros((D, D), np.complex128)
    for a in range(D):
        for b in range(D):
            if abs(h[b, a]) > tol:
                outs[a, cnt[a]] = b
                vals[a, cnt[a]] = h[b, a]
                cnt[a] += 1
    return cnt, outs, vals


@njit(cache=True)
def _build_rows(L, r, cnt, outs, vals, reps, per, oid, shift, src, inv, m):
    """Row c of M = H^T (M[c, row] = <row|H|c>) for the momentum sector m.
    Diagonal contributions are accumulated into a single entry per row."""
    n = len(src)
    maxper = 0
    for a in range(cnt.shape[0]):
        if cnt[a] > maxper:
            maxper = cnt[a]
    cap = n * (L * maxper + 1)
    indptr = np.zeros(n + 1, np.int64)
    indices = np.empty(cap, np.int32)
    data = np.empty(cap, np.complex128)
    k = 0
    twopi = 2.0 * np.pi
    for c in range(n):
        rr = src[c]
        s = reps[rr]
        diag = 0.0 + 0.0j
        for j in range(L):
            a = 0
            for i in range(r):
                a = (a << 1) | ((s >> ((j + i) % L)) & 1)
            for q in range(cnt[a]):
                b = outs[a, q]
                if b == a:
                    diag += vals[a, q]
                    continue
                sp = s
                for i in range(r):
                    bit_new = (b >> (r - 1 - i)) & 1
                    pos = (j + i) % L
                    if ((sp >> pos) & 1) != bit_new:
                        sp ^= (1 << pos)
                o = oid[sp]
                row = inv[o]
                if row < 0:
                    continue
                ph = twopi * m * shift[sp] / L
                indices[k] = row
                # store conj(<row|H|c>) so that the CSR rows directly give H (Hermitian)
                z = vals[a, q] * np.sqrt(per[rr] / per[o]) * (np.cos(ph) + 1j * np.sin(ph))
                data[k] = z.conjugate()
                k += 1
        indices[k] = c
        data[k] = diag.conjugate()
        k += 1
        indptr[c + 1] = k
    return indptr, indices[:k], data[:k]


def sector_hamiltonian(L, h, m):
    r = int(round(np.log2(h.shape[0])))
    reps, per, oid, shift = orbits(L)
    src, inv = sector(L, m)
    cnt, outs, vals = _local_table(h)
    indptr, indices, data = _build_rows(L, r, cnt, outs, vals, reps, per, oid, shift, src, inv, m)
    if np.max(np.abs(data.imag), initial=0) < 1e-13:
        data = data.real.copy()
    # row c of this matrix is conj(column c of H) = row c of H, since H is Hermitian
    return sparse.csr_matrix((data, indices, indptr), shape=(len(src), len(src)), copy=False)


def lowest(H, k=4, tol=1e-10, seed=7):
    n = H.shape[0]
    if n <= 400:
        e, u = np.linalg.eigh(H.toarray())
        return e[:k], u[:, :k]
    rng = np.random.default_rng(seed)
    v0 = rng.normal(size=n)
    if np.iscomplexobj(H.data):
        v0 = v0 + 1j * rng.normal(size=n)
    k = min(k, n - 2)
    e, u = eigsh(H, k=k, which='SA', tol=tol, v0=v0, ncv=max(3 * k + 4, 30), maxiter=50000)
    ix = np.argsort(e)
    return e[ix], u[:, ix]


def expand(L, m, u):
    """Sector vector -> full 2^L vector (normalized)."""
    reps, per, oid, shift = orbits(L)
    src, inv = sector(L, m)
    idx = inv[oid]
    good = idx >= 0
    psi = np.zeros(1 << L, complex)
    psi[good] = u[idx[good]] * np.exp(-2j * np.pi * m * shift[good] / L) / np.sqrt(per[oid[good]].astype(float))
    return psi / np.linalg.norm(psi)


def entropies(psi, L, lmax=None):
    """Von Neumann entropy (natural log) of the block of the first l sites, l=1..lmax."""
    lmax = L // 2 if lmax is None else lmax
    out = []
    for l in range(1, lmax + 1):
        M = psi.reshape(1 << (L - l), 1 << l)
        rho = M.conj().T @ M
        w = np.linalg.eigvalsh(rho)
        w = w[w > 1e-15]
        out.append(float(-(w * np.log(w)).sum()))
    return np.array(out)


def full_hamiltonian(L, h):
    r = int(round(np.log2(h.shape[0])))
    Hd = np.zeros((1 << L, 1 << L), complex)
    eye = np.eye(2)
    for j in range(L):
        # build operator on sites j..j+r-1 by permutation of a dense kron
        ops = h
        # embed: index ordering with site L-1 most significant in kron
        full = np.zeros((1 << L, 1 << L), complex)
        for s in range(1 << L):
            a = 0
            for i in range(r):
                a = (a << 1) | ((s >> ((j + i) % L)) & 1)
            for b in range(1 << r):
                v = ops[b, a]
                if v == 0:
                    continue
                sp = s
                for i in range(r):
                    pos = (j + i) % L
                    nb = (b >> (r - 1 - i)) & 1
                    if ((sp >> pos) & 1) != nb:
                        sp ^= 1 << pos
                full[sp, s] += v
        Hd += full
    return Hd


def spectrum(L, h, k=4, sectors=None, want_vec=False):
    """Lowest k levels per momentum sector m (k = 2 pi m / L)."""
    sectors = range(L) if sectors is None else sectors
    out = {}
    vecs = {}
    for m in sectors:
        H = sector_hamiltonian(L, h, m)
        e, u = lowest(H, k=k)
        out[m] = e
        if want_vec:
            vecs[m] = u[:, 0]
    return (out, vecs) if want_vec else out


def selftest():
    rng = np.random.default_rng(1)
    ok = True
    for L, r in ((6, 2), (7, 2), (6, 3), (8, 2)):
        A = rng.normal(size=(1 << r, 1 << r)) + 1j * rng.normal(size=(1 << r, 1 << r))
        h = A + A.conj().T
        Ef = np.sort(np.linalg.eigvalsh(full_hamiltonian(L, h)))
        Es = []
        for m in range(L):
            H = sector_hamiltonian(L, h, m).toarray()
            assert np.allclose(H, H.conj().T)
            Es.extend(np.linalg.eigvalsh(H))
        Es = np.sort(Es)
        err = np.max(np.abs(Es - Ef))
        print(f"selftest L={L} r={r}: max spectrum error {err:.2e}")
        ok &= err < 1e-10
        # momentum eigenstate check: expand and verify H psi = E psi in full space
        Hf = full_hamiltonian(L, h)
        for m in (0, 1, L - 1):
            e, u = np.linalg.eigh(sector_hamiltonian(L, h, m).toarray())
            psi = expand(L, m, u[:, 0])
            res = np.linalg.norm(Hf @ psi - e[0] * psi)
            ok &= res < 1e-9
            # translation eigenvalue
            s = np.arange(1 << L)
            dest = ((s << 1) & ((1 << L) - 1)) | (s >> (L - 1))
            Tpsi = np.empty_like(psi); Tpsi[dest] = psi
            ok &= abs(np.vdot(psi, Tpsi) - np.exp(2j * np.pi * m / L)) < 1e-9 or abs(abs(np.vdot(psi, Tpsi)) - 1) < 1e-9
    print("selftest", "PASSED" if ok else "FAILED")
    return ok


if __name__ == '__main__':
    selftest()

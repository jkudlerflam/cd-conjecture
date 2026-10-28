#!/usr/bin/env python3
"""Momentum- and parity-resolved exact diagonalization of qubit rings.

H = sum_j T^j h T^-j, h an r-site term (r = 2 here, r = 3 for the cluster-Ising
control).  Site j is bit j of the basis label, bit value 1 means sigma^z = -1,
and inside h the leftmost site is the most significant kron factor.  Sectors
are labelled by the momentum k = 2 pi m / L and, when the bond commutes with
P = prod_j sigma^z_j, by P = +1 or -1 (even or odd number of 1 bits).

The orbit construction and the CSR builder are copied from
../d2_classification/ed.py (Claude Opus 5.5, 2026-10-27).  New here: the parity
filter, the model parametrizations, and the helpers used by the scans.
selftest() compares every (m, P) sector with dense full-space diagonalization.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
from functools import lru_cache
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from numba import njit

X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
PAULI = np.array([X, Y, Z])
I2 = np.eye(2)


# --------------------------------------------------------------------------- bonds
def bond(J, field):
    """h = sum_ab J_ab s^a (x) s^b + sum_a f_a (s^a (x) 1 + 1 (x) s^a)/2, so that
    H = sum_j [s_j^T J s_{j+1} + f . s_j]  (same convention as d2_classification)."""
    J = np.asarray(J, float); field = np.asarray(field, float)
    h = np.zeros((4, 4), complex)
    for a in range(3):
        h += field[a] * 0.5 * (np.kron(PAULI[a], I2) + np.kron(I2, PAULI[a]))
        for b in range(3):
            h += J[a, b] * np.kron(PAULI[a], PAULI[b])
    return h


def z2_params(Jx, Jy, Jz, D=0.0, h=0.0, hx=0.0, Jxz=0.0):
    """C15 convention H = sum [Jx XX + Jy YY + Jz ZZ + D (XY - YX) - h Z]
    plus the Z2-breaking extensions hx X (field) and Jxz (XZ + ZX).
    Returns (J, field) in the bond() convention."""
    J = np.diag([Jx, Jy, Jz]).astype(float)
    J[0, 1] += D; J[1, 0] -= D
    J[0, 2] += Jxz; J[2, 0] += Jxz
    return J, np.array([hx, 0.0, -h])


def kitaev(tau=1.0, delta=0.0, V=0.0, mu=0.0, D=0.0, hx=0.0, Jxz=0.0):
    """Interacting Kitaev chain (C15 eq. 1, real tau) plus a DM/complex-hopping D:
    tau = -(Jx+Jy), delta = Jx-Jy, V = 4 Jz, mu = V - 2h."""
    Jx = (-tau + delta) / 2; Jy = (-tau - delta) / 2; Jz = V / 4; h = (V - mu) / 2
    return z2_params(Jx, Jy, Jz, D, h, hx, Jxz)


# --------------------------------------------------------------------------- orbits
@njit(cache=True)
def _orbits(L):
    N = 1 << L
    oid = np.full(N, -1, np.int32)
    shift = np.zeros(N, np.int8)
    reps = np.empty(N, np.int64)
    per = np.empty(N, np.int16)
    par = np.empty(N, np.int8)
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
        c = 0
        x = s
        while x:
            c += x & 1
            x >>= 1
        par[nr] = 1 - 2 * (c & 1)
        nr += 1
    return reps[:nr].copy(), per[:nr].copy(), oid, shift, par[:nr].copy()


@lru_cache(maxsize=3)
def orbits(L):
    return _orbits(L)


def sector(L, m, P=None):
    reps, per, oid, shift, par = orbits(L)
    allowed = (per.astype(np.int64) * m) % L == 0
    if P is not None:
        allowed &= par == P
    active = np.flatnonzero(allowed).astype(np.int64)
    inv = np.full(len(reps), -1, np.int64)
    inv[active] = np.arange(len(active))
    return active, inv


def _local_table(h, tol=1e-14):
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
                z = vals[a, q] * np.sqrt(per[rr] / per[o]) * (np.cos(ph) + 1j * np.sin(ph))
                data[k] = z.conjugate()
                k += 1
        indices[k] = c
        data[k] = diag.conjugate()
        k += 1
        indptr[c + 1] = k
    return indptr, indices[:k], data[:k]


def has_parity(h, tol=1e-12):
    """True if the bond commutes with Z (x) Z (so H commutes with prod Z)."""
    r = int(round(np.log2(h.shape[0])))
    P = Z
    for _ in range(r - 1):
        P = np.kron(P, Z)
    return np.max(np.abs(P @ h - h @ P)) < tol


def sector_hamiltonian(L, h, m, P=None):
    r = int(round(np.log2(h.shape[0])))
    reps, per, oid, shift, par = orbits(L)
    src, inv = sector(L, m, P)
    cnt, outs, vals = _local_table(h)
    indptr, indices, data = _build_rows(L, r, cnt, outs, vals, reps, per, oid, shift, src, inv, m)
    if np.max(np.abs(data.imag), initial=0) < 1e-13:
        data = data.real.copy()
    return sparse.csr_matrix((data, indices, indptr), shape=(len(src), len(src)), copy=False)


def lowest(H, k=4, tol=1e-10, seed=7):
    n = H.shape[0]
    if n == 0:
        return np.full(k, np.nan), np.zeros((0, k))
    if n <= 600:
        e, u = np.linalg.eigh(H.toarray())
        kk = min(k, n)
        e2 = np.full(k, np.nan); e2[:kk] = e[:kk]
        return e2, u[:, :kk]
    rng = np.random.default_rng(seed)
    v0 = rng.normal(size=n)
    if np.iscomplexobj(H.data):
        v0 = v0 + 1j * rng.normal(size=n)
    kk = min(k, n - 2)
    e, u = eigsh(H, k=kk, which='SA', tol=tol, v0=v0, ncv=max(3 * kk + 6, 30), maxiter=50000)
    ix = np.argsort(e)
    return e[ix], u[:, ix]


def expand(L, m, u, P=None):
    reps, per, oid, shift, par = orbits(L)
    src, inv = sector(L, m, P)
    idx = inv[oid]
    good = idx >= 0
    psi = np.zeros(1 << L, complex)
    psi[good] = u[idx[good]] * np.exp(-2j * np.pi * m * shift[good] / L) / np.sqrt(per[oid[good]].astype(float))
    return psi / np.linalg.norm(psi)


def entropies(psi, L, lmax=None):
    lmax = L // 2 if lmax is None else lmax
    out = []
    for l in range(1, lmax + 1):
        M = psi.reshape(1 << (L - l), 1 << l)
        s = np.linalg.svd(M, compute_uv=False) ** 2
        s = s[s > 1e-15]
        out.append(float(-(s * np.log(s)).sum()))
    return np.array(out)


def c_chord(S, L, lmin=None):
    """PBC fit S(l) = (c/3) ln[(L/pi) sin(pi l/L)] + b over l >= L/4."""
    lmin = int(np.ceil(L / 4)) if lmin is None else lmin
    ls = np.arange(lmin, L // 2 + 1)
    x = np.log(L / np.pi * np.sin(np.pi * ls / L)) / 3
    return float(np.polyfit(x, np.asarray(S)[ls - 1], 1)[0])


def spectrum(L, h, k=3, parity=None, sectors=None, want_gs=False):
    """Lowest k levels in every (m, P) sector.  Returns dict {(m,P): energies};
    P is 0 when the parity is not resolved.  With want_gs, also the ground state
    sector and its full vector."""
    if parity is None:
        parity = has_parity(h)
    Ps = (1, -1) if parity else (None,)
    sectors = range(L) if sectors is None else sectors
    out, gs = {}, None
    for m in sectors:
        for P in Ps:
            H = sector_hamiltonian(L, h, m, P)
            e, u = lowest(H, k=k)
            key = (m, P if P is not None else 0)
            out[key] = e
            if want_gs and (gs is None or e[0] < gs[0] - 1e-12):
                gs = (e[0], key, u[:, 0])
    if want_gs:
        e0, key, u = gs
        psi = expand(L, key[0], u, key[1] if key[1] != 0 else None)
        return out, key, psi
    return out


def full_dense(L, h):
    r = int(round(np.log2(h.shape[0])))
    Hd = np.zeros((1 << L, 1 << L), complex)
    for j in range(L):
        for s in range(1 << L):
            a = 0
            for i in range(r):
                a = (a << 1) | ((s >> ((j + i) % L)) & 1)
            for b in range(1 << r):
                v = h[b, a]
                if v == 0:
                    continue
                sp = s
                for i in range(r):
                    pos = (j + i) % L
                    nb = (b >> (r - 1 - i)) & 1
                    if ((sp >> pos) & 1) != nb:
                        sp ^= 1 << pos
                Hd[sp, s] += v
    return Hd


def selftest():
    rng = np.random.default_rng(3)
    ok = True
    for L in (6, 7, 8):
        # a Z2-symmetric random bond (C15 family) and a generic one
        for hb in (bond(*z2_params(*rng.normal(size=5))), None):
            if hb is None:
                A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)); hb = A + A.conj().T
            Ef = np.sort(np.linalg.eigvalsh(full_dense(L, hb)))
            par = has_parity(hb)
            Es = []
            for m in range(L):
                for P in ((1, -1) if par else (None,)):
                    Es.extend(np.linalg.eigvalsh(sector_hamiltonian(L, hb, m, P).toarray()))
            err = np.max(np.abs(np.sort(Es) - Ef))
            print(f'selftest L={L} parity={par}: max error {err:.1e}')
            ok &= err < 1e-10
    print('selftest', 'PASSED' if ok else 'FAILED')
    return ok


if __name__ == '__main__':
    selftest()

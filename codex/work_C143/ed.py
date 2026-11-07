"""Momentum-resolved exact diagonalization for a uniform nearest-neighbour qubit ring.

Bond convention (as in specs/_d2_background.md):
    h = sum_{ab} J[a,b] s^a (x) s^b + sum_a B[a] (s^a (x) 1 + 1 (x) s^a)/2,
with Pauli matrices s^a, a = x, y, z. Ring H_L = sum_i h_{i,i+1}, periodic.
Basis bit convention: site i is bit (L-1-i) of the integer (site 0 is the most significant bit).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [SX, SY, SZ]


def bond(J, B=(0, 0, 0), stag=(0, 0, 0)):
    """4x4 bond matrix. stag adds the gauge term o.(s(x)1 - 1(x)s), which cancels on rings."""
    J = np.asarray(J, float).reshape(3, 3)
    h = np.zeros((4, 4), complex)
    for a in range(3):
        for b in range(3):
            h += J[a, b] * np.kron(PAULI[a], PAULI[b])
        h += B[a] * (np.kron(PAULI[a], I2) + np.kron(I2, PAULI[a])) / 2
        h += stag[a] * (np.kron(PAULI[a], I2) - np.kron(I2, PAULI[a]))
    return h


def params_to_JB(x):
    x = np.asarray(x, float)
    return x[:9].reshape(3, 3), x[9:12]


def gauge_M(h):
    """M = min over gauge o of (lmax - lmin)/2 of h + o.(s(x)1 - 1(x)s)."""
    G = [np.kron(P, I2) - np.kron(I2, P) for P in PAULI]

    def f(o):
        w = np.linalg.eigvalsh(h + sum(o[a] * G[a] for a in range(3)))
        return (w[-1] - w[0]) / 2

    best = None
    for o0 in [np.zeros(3), 0.3 * np.ones(3), -0.3 * np.ones(3)]:
        r = minimize(f, o0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-9, maxiter=2000))
        if best is None or r.fun < best.fun:
            best = r
    return best.fun, best.x


def gauge_lmin(h):
    """max over gauge o of lambda_min(h + gauge): the frustration-free reference energy per bond."""
    G = [np.kron(P, I2) - np.kron(I2, P) for P in PAULI]

    def f(o):
        return -np.linalg.eigvalsh(h + sum(o[a] * G[a] for a in range(3)))[0]

    best = None
    for o0 in [np.zeros(3), 0.3 * np.ones(3), -0.3 * np.ones(3), np.array([0.5, -0.2, 0.1])]:
        r = minimize(f, o0, method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-12, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
    return -best.fun, best.x


_cache = {}


def _orbits(L):
    if L in _cache:
        return _cache[L]
    N = 1 << L
    x = np.arange(N, dtype=np.int64)
    mask = N - 1
    rots = np.empty((L, N), dtype=np.int64)
    cur = x.copy()
    for j in range(L):
        rots[j] = cur
        # translate by one site: site i -> i+1, i.e. rotate bits right by one
        cur = ((cur >> 1) | ((cur & 1) << (L - 1))) & mask
    rep = rots.min(axis=0)
    # shift l such that x = T^l rep  <=> rep = T^{-l} x ; find j with rots[j] == rep, then rep = T^j x, so x = T^{-j} rep
    j = rots.argmin(axis=0)
    shift = (-j) % L
    # period of each orbit
    per = np.full(N, L, dtype=np.int64)
    for p in range(1, L):
        if L % p == 0:
            hit = (rots[p] == x) & (per == L)
            per[hit] = p
    _cache[L] = (rep, shift, per)
    return _cache[L]


def ham_sparse(h, L):
    N = 1 << L
    x = np.arange(N, dtype=np.int64)
    rows, cols, vals = [], [], []
    for i in range(L):
        a = L - 1 - i
        b = L - 1 - ((i + 1) % L)
        bi = (x >> a) & 1
        bj = (x >> b) & 1
        inidx = 2 * bi + bj
        base = x & ~((1 << a) | (1 << b))
        for out in range(4):
            amp = h[out, inidx]
            nz = np.abs(amp) > 1e-14
            if not nz.any():
                continue
            oi, oj = out >> 1, out & 1
            y = base | (oi << a) | (oj << b)
            rows.append(y[nz]); cols.append(x[nz]); vals.append(amp[nz])
    rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N))


def momentum_projector(L, m):
    rep, shift, per = _orbits(L)
    N = 1 << L
    k = 2 * np.pi * m / L
    ok = ((m * per) % L) == 0
    reps = np.unique(rep[ok])
    col = -np.ones(N, dtype=np.int64)
    idx = np.searchsorted(reps, rep[ok])
    col[ok] = idx
    xs = np.nonzero(ok)[0]
    vals = np.exp(-1j * k * shift[xs]) / np.sqrt(per[xs])
    P = sp.csr_matrix((vals, (xs, col[xs])), shape=(N, len(reps)))
    return P


def spectrum(h, L, nev=None, dense_max=2500):
    """Return dict m -> sorted eigenvalues (all if dense, else nev lowest)."""
    H = ham_sparse(h, L)
    out = {}
    for m in range(L):
        P = momentum_projector(L, m)
        if P.shape[1] == 0:
            out[m] = np.array([])
            continue
        Hk = (P.conj().T @ H @ P)
        n = Hk.shape[0]
        if n <= dense_max or nev is None or nev >= n - 2:
            w = np.linalg.eigvalsh(Hk.toarray())
            out[m] = w if nev is None else w[:nev]
        else:
            w = spla.eigsh(Hk, k=nev, which="SA", tol=1e-10, return_eigenvectors=False)
            out[m] = np.sort(w)
    return out


def all_levels(spec):
    ev, ks = [], []
    for m, w in spec.items():
        ev.append(w); ks.append(np.full(len(w), m))
    return np.concatenate(ev), np.concatenate(ks)


if __name__ == "__main__":
    # self-tests against exact results
    # 1) XX chain h = (XX+YY)/2 ... compare with free fermions on L=8 (both parities)
    L = 8
    h = bond(np.diag([1, 1, 0]))
    w, _ = all_levels(spectrum(h, L))
    Hd = ham_sparse(h, L).toarray()
    wd = np.linalg.eigvalsh(Hd)
    print("XX L=8 momentum vs dense max diff", np.max(np.abs(np.sort(w) - wd)))
    # 2) Heisenberg h = XX+YY+ZZ, E0(L=8) known = -3.651093408937 * 4 (for S.S): E0(S.S, L=8) = -3.651093
    h = bond(np.eye(3))
    w, _ = all_levels(spectrum(h, 8))
    print("Heisenberg L=8 E0/4 =", w.min() / 4, "(expect -3.6510934)")
    # 3) TFIM critical: h = -XX - (Z1+Z2)/2 -> H = -sum XX - sum Z ; E0(L) = -sum_k 2|cos(k/2)| over antiperiodic k
    L = 10
    h = bond(-np.diag([1, 0, 0]), B=(0, 0, -1))
    w, _ = all_levels(spectrum(h, L))
    ks = np.pi * (2 * np.arange(L) + 1) / L
    print("TFIM L=10 E0", w.min(), "exact", -np.sum(2 * np.abs(np.cos(ks / 2))))
    print("M(Heis)", gauge_M(bond(np.eye(3)))[0], "(expect 2 for S=sigma normalization: (1-(-3))/2)")

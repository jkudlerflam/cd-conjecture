"""Parity-symmetric qubit class: bond h = Jx XX + Jy YY + D (XY - YX) + Jz ZZ + (hz/2)(Z1 + Z2).

Ring and open-chain ED in parity sectors (Pi = prod Z), entanglement, free-fermion (V = Jz = 0) BdG dispersion.
Conventions: |0> = spin up (Z = +1). Jordan-Wigner with Z_j = 1 - 2 n_j gives
  H = sum_j [ w c+_{j+1} c_j + h.c. - (Jx - Jy)(c_j c_{j+1} + h.c.) + 4 Jz n_j n_{j+1} - (4 Jz + 2 hz) n_j ] + const,
  w = Jx + Jy + 2 i D.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def bond(Jx, Jy, D, Jz, hz):
    h = (Jx * np.kron(SX, SX) + Jy * np.kron(SY, SY) + D * (np.kron(SX, SY) - np.kron(SY, SX))
         + Jz * np.kron(SZ, SZ) + 0.5 * hz * (np.kron(SZ, I2) + np.kron(I2, SZ)))
    return h


def is_real(h):
    return np.max(np.abs(h.imag)) < 1e-14


def sector_states(L, parity):
    s = np.arange(2 ** L, dtype=np.int64)
    pc = np.zeros_like(s)
    t = s.copy()
    for _ in range(L):
        pc += t & 1; t >>= 1
    sel = s[(pc & 1) == parity]
    idx = -np.ones(2 ** L, dtype=np.int64)
    idx[sel] = np.arange(len(sel))
    return sel, idx


def ham(h, L, parity, periodic=True, sector=None):
    """Sparse H in the Pi-sector with parity of the number of down spins = parity. Site j <-> bit (L-1-j)."""
    if sector is None:
        sector = sector_states(L, parity)
    sel, idx = sector
    n = len(sel)
    h = np.asarray(h)
    dt = float if is_real(h) else complex
    hh = h.real if dt is float else h
    rows, cols, vals = [], [], []
    nb = L if periodic else L - 1
    for j in range(nb):
        a, b = j, (j + 1) % L
        ba, bb = L - 1 - a, L - 1 - b
        sa = (sel >> ba) & 1; sb = (sel >> bb) & 1
        loc = 2 * sa + sb
        for lp in range(4):
            for l in range(4):
                v = hh[lp, l]
                if abs(v) < 1e-15:
                    continue
                m = loc == l
                if not m.any():
                    continue
                st = sel[m]
                na, nb_ = lp >> 1, lp & 1
                new = (st & ~((1 << ba) | (1 << bb))) | (na << ba) | (nb_ << bb)
                rows.append(idx[new]); cols.append(np.nonzero(m)[0]); vals.append(np.full(m.sum(), v, dtype=dt))
    rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    return H


def lowest(H, k=2):
    n = H.shape[0]
    if n <= 600:
        w, v = np.linalg.eigh(H.toarray())
        return w[:k], v[:, :k]
    w, v = spla.eigsh(H, k=k, which="SA", tol=1e-10, ncv=max(20, 2 * k + 5))
    o = np.argsort(w)
    return w[o], v[:, o]


def half_entropy(vec, sel, L, l=None):
    if l is None:
        l = L // 2
    full = np.zeros(2 ** L, dtype=vec.dtype)
    full[sel] = vec
    M = full.reshape(2 ** l, 2 ** (L - l))
    s = np.linalg.svd(M, compute_uv=False)
    p = s ** 2; p = p[p > 1e-15]; p /= p.sum()
    return float(-np.sum(p * np.log(p)))


def ring_data(h, L, k=3):
    """Lowest k levels in each parity sector, global ground state's half-ring entropy."""
    out = {}
    best = None
    for par in (0, 1):
        sec = sector_states(L, par)
        H = ham(h, L, par, True, sec)
        w, v = lowest(H, k)
        out[par] = w
        if best is None or w[0] < best[0]:
            best = (w[0], v[:, 0], sec[0])
    S = half_entropy(best[1], best[2], L)
    lev = np.sort(np.concatenate([out[0], out[1]]))
    return dict(E0=float(lev[0]), levels=[float(x) for x in lev], S=S, E0_par=[float(out[0][0]), float(out[1][0])])


def free_f(k, Jx, Jy, D, hz):
    """Free (Jz = 0) Bogoliubov energy E_+(k) of one quasiparticle band (C143 eq. 2.1 with a = 2D,
    c = Jx + Jy, g = |Jx - Jy|, B = -hz up to sign conventions): E = 2 f."""
    a = 2 * D; c = Jx + Jy; g = abs(Jx - Jy); B = hz
    return 2 * (-a * np.sin(k) + np.sqrt((c * np.cos(k) + B) ** 2 + g ** 2 * np.sin(k) ** 2))


def spin_from_fermion(t, phi, Delta, V, mu):
    """Fermion form sum[w c+_{j+1}c_j + h.c. - Delta(c_j c_{j+1} + h.c.) + V n_j n_{j+1} - mu n_j], w = t e^{i phi}
    -> spin couplings (Jx, Jy, D, Jz, hz). Free (V = 0) Bogoliubov band:
    E_+(k) = 2 t sin(phi) sin k + sqrt((2 t cos(phi) cos k - mu)^2 + 4 Delta^2 sin^2 k)."""
    c = t * np.cos(phi); a = t * np.sin(phi)
    Jx = (c + Delta) / 2; Jy = (c - Delta) / 2; D = a / 2; Jz = V / 4; hz = (mu - V) / 2
    return (Jx, Jy, D, Jz, hz)


def free_band(k, t, phi, Delta, mu):
    return 2 * t * np.sin(phi) * np.sin(k) + np.sqrt((2 * t * np.cos(phi) * np.cos(k) - mu) ** 2 + 4 * Delta ** 2 * np.sin(k) ** 2)

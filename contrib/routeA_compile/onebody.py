"""Exact second-order one-body logical Hamiltonian of a two-cell junction,
and the Fermi-point count of the resulting n-band logical chain.

Total fermion number is exact, so the order-lambda^2 logical bond restricted
to the 0- and 1-data-fermion code sectors is E0 + sum h_ab c_a^dag c_b
exactly (whatever the higher-body parts are).  Intermediate states from a
<=1-data-fermion code state carry at most 3 data fermions and 2 carriers,
so the truncation K=2, max_data=3 is exact for h.  H0 conserves the carrier
number of each cell and V puts exactly one carrier in each cell, so the
basis {code} + {one carrier per cell} is also exact for h.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
import scipy.sparse.linalg as sla
from model import Chain, parts, assemble, DEFAULT

_CACHE = {}


def prepared(n):
    if n not in _CACHE:
        ch = Chain(n, 2, K=2, max_data=3, one_per_cell=True)
        _CACHE[n] = (ch, parts(ch))
    return _CACHE[n]


def logical_h(n, params=None):
    ch, P = prepared(n)
    p = dict(DEFAULT)
    if params:
        p.update(params)
    H0, V = assemble(P, p)
    dm = ch.data_modes
    # code states with <=1 data fermion
    vac = [c for c in ch.code if all(((ch.states[c] >> m) & 1) == 0 for m in dm)][0]
    kets = []   # (index, sign) of c_m^dag |vac>
    for m in dm:
        sg, new = Chain.apply([(m, 1)], [], ch.states[vac])
        kets.append((ch.index[new], sg))
    cols = [vac] + [k for k, _ in kets]
    signs = np.array([1] + [s for _, s in kets], float)
    mask = np.ones(ch.dim, bool)
    mask[ch.code] = False
    Q = np.where(mask)[0]
    VQP = V[Q][:, cols].toarray() * signs[None, :]
    X = sla.spsolve(H0[Q][:, Q].tocsc(), VQP)
    H2 = -VQP.conj().T @ X
    H2 = 0.5 * (H2 + H2.conj().T)
    E0 = H2[0, 0].real
    h = H2[1:, 1:] - E0 * np.eye(len(dm))
    # zero-to-one-fermion blocks must vanish by number conservation
    assert np.max(abs(H2[0, 1:])) < 1e-14
    return E0, h


def blocks(h, n):
    T = h[:n, n:]                     # c_{j,a}^dag T_ab c_{j+1,b}
    M = h[:n, :n] + h[n:, n:]         # each cell has a right and a left junction
    return T, M


def bands(T, M, nk=4000):
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    return np.array([np.linalg.eigvalsh(M + T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k)) for k in ks])


def fermi_count(E, mu):
    """Number of Fermi points at chemical potential mu (sorted bands E[k, i])."""
    s = E < mu
    return int(np.sum(s != np.roll(s, 1, axis=0)))


def max_c(T, M, nk=4000):
    """Largest c = (#Fermi points)/2 over all chemical potentials.

    c(mu) only changes at band extrema, so test the midpoints between all
    consecutive extremal values of the sorted bands on a fine k grid."""
    scale = max(np.abs(T).max(), np.abs(M).max())
    T, M = T / scale, M / scale          # c(mu) is scale invariant
    E = bands(T, M, nk)
    ext = []
    for i in range(E.shape[1]):
        e = E[:, i]
        lo = (e <= np.roll(e, 1)) & (e <= np.roll(e, -1))
        hi = (e >= np.roll(e, 1)) & (e >= np.roll(e, -1))
        ext += list(e[lo | hi])
    ext = np.unique(np.round(np.array(ext), 15))
    best, arg = 0, None
    for mu in 0.5 * (ext[1:] + ext[:-1]):
        f = fermi_count(E, mu)
        if f > best:
            best, arg = f, float(mu) * scale
    return best / 2, arg

"""Charges and product-state configurations for qdmrg.py runs.  Claude (Opus 5.5), 2026-10-28."""
import numpy as np


def lv_charges(n):
    """Local states of LV_n: 0 = vacancy, 1 + b (bit f of b, flavour 0 most significant) = cell.
    Charge vector (M, N_0, ..., N_{n-1})."""
    q = [[0] * (n + 1)]
    for b in range(2 ** n):
        q.append([1] + [(b >> (n - 1 - f)) & 1 for f in range(n)])
    return np.array(q)


def spread(L, N):
    """Boolean array of N evenly spread ones on L sites."""
    x = np.zeros(L, bool)
    if N:
        x[np.floor((np.arange(N) + 0.5) * L / N).astype(int)] = True
    return x


def lv_config(L, M, N):
    n = len(N)
    occ = spread(L, M)
    fl = [spread(M, Nf) for Nf in N]
    if n > 1:   # stagger the second flavour so cells are not doubly occupied first
        fl[1] = np.roll(fl[1], 1)
    cfg = []
    m = 0
    for j in range(L):
        if not occ[j]:
            cfg.append(0); continue
        b = 0
        for f in range(n):
            b |= int(fl[f][m]) << (n - 1 - f)
        cfg.append(1 + b); m += 1
    return cfg


def kxx_bond(k, vels=None):
    """k decoupled XX chains with hoppings vels (default 1) on a d = 2^k site (state bit i = occupation of
    chain i, chain 0 most significant).  The exact ground-state entanglement of a fixed-N sector does not
    depend on the hoppings; a truncated (finite-chi) DMRG state does."""
    vels = [1.0] * k if vels is None else list(vels)
    sp = np.array([[0., 1.], [0., 0.]])      # |0><1| annihilates
    I = np.eye(2)

    def emb(o, i):
        m = np.ones((1, 1))
        for j in range(k):
            m = np.kron(m, o if j == i else I)
        return m
    h = 0
    for i in range(k):
        A = np.kron(emb(sp.T, i), emb(sp, i))
        h = h - vels[i] * (A + A.T)
    q = np.array([[(b >> (k - 1 - i)) & 1 for i in range(k)] for b in range(2 ** k)])
    return h, q


def kxx_config(L, N):
    k = len(N)
    fl = [np.roll(spread(L, Nf), i) for i, Nf in enumerate(N)]
    return [sum(int(fl[i][j]) << (k - 1 - i) for i in range(k)) for j in range(L)]


def ff_open_S(L, N):
    """Exact S(l), l = 0..L, of the open XX chain (free fermions) ground state with N particles."""
    if N == 0 or N == L:
        return np.zeros(L + 1)
    x = np.arange(1, L + 1)
    k = np.arange(1, N + 1)
    phi = np.sqrt(2.0 / (L + 1)) * np.sin(np.pi * np.outer(x, k) / (L + 1))
    C = phi @ phi.T
    out = [0.0]
    for l in range(1, L):
        nu = np.linalg.eigvalsh(C[:l, :l]).clip(1e-15, 1 - 1e-15)
        out.append(float(-(nu * np.log(nu) + (1 - nu) * np.log(1 - nu)).sum()))
    return np.array(out + [0.0])


def ref_S(L, Ns):
    """Exact reference: sum of independent open free-fermion chains with particle numbers Ns."""
    return sum(ff_open_S(L, N) for N in Ns)

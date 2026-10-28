"""Independent wrapper pieces: vacancy charge fluid, endpoint helices,
compensator, selector head, and the wrapped two-cell second-order bond.

Written from codex/work_C44/PROOFS.md Sections 1, 4, 5 and the C8 notes
length_filter.md, arrow_wrapper.md, selector_wrapper.md, helix_equalizer.md.
"""
import itertools
import numpy as np
import scipy.sparse as sp
from core46 import Core, Params, SPECIES, port_terms


# ---------------------------------------------------------------- charge fluid
class Charge:
    """k occupied labels on N=k+p physical sites, ordered coordinates."""

    def __init__(self, k, p, U, wvac):
        self.k, self.p, self.N = k, p, k + p
        N = self.N
        self.configs = [c for c in itertools.combinations(range(N), k)]
        self.index = {c: i for i, c in enumerate(self.configs)}
        m = len(self.configs)
        H = np.zeros((m, m))
        for i, c in enumerate(self.configs):
            occ = set(c)
            H[i, i] = U * sum(1 for x in range(N - 1) if (x in occ) != (x + 1 in occ))
            for j, o in enumerate(c):
                for o2 in (o - 1, o + 1):
                    if 0 <= o2 < N and o2 not in occ:
                        c2 = tuple(sorted((occ - {o}) | {o2}))
                        H[self.index[c2], i] += -wvac
        self.H = H
        # diagonalize in reflection-even and -odd sectors separately: the
        # tunnelling doublet is split by ~U rho^(kp) and a plain eigh mixes it
        refl = [self.index[tuple(sorted(N - 1 - x for x in c))] for c in self.configs]
        ev, od, seen = [], [], set()
        for i, j in enumerate(refl):
            if i in seen:
                continue
            seen |= {i, j}
            v = np.zeros(m)
            if i == j:
                v[i] = 1.0
                ev.append(v)
            else:
                v[i] = v[j] = 1 / np.sqrt(2)
                ev.append(v)
                w = np.zeros(m)
                w[i], w[j] = 1 / np.sqrt(2), -1 / np.sqrt(2)
                od.append(w)
        Be, Bo = np.array(ev).T, np.array(od).T
        Ee, Ve = np.linalg.eigh(Be.T @ H @ Be)
        Eo, Vo = np.linalg.eigh(Bo.T @ H @ Bo)
        E = np.concatenate([Ee, Eo])
        V = np.hstack([Be @ Ve, Bo @ Vo])
        order = np.argsort(E, kind='stable')
        E, V = E[order], V[:, order]
        self.eig = (E, V)
        self.e = E[0]
        self.gap = E[1] - E[0]
        chi = V[:, 0] * np.sign(V[:, 0].sum())
        assert (chi > 0).all() and abs(Ee[0] - E[0]) == 0
        self.chi = chi
        self.reflection_asymmetry = float(np.abs(chi - chi[refl]).max())
        self.Dj = np.array([[1.0 if c[j + 1] == c[j] + 1 else 0.0 for j in range(k - 1)]
                            for c in self.configs])
        self.OL = np.array([1.0 if c[0] == 0 else 0.0 for c in self.configs])
        self.OR = np.array([1.0 if c[-1] == N - 1 else 0.0 for c in self.configs])
        self.r = float(chi @ (self.OL * chi))
        self.rR = float(chi @ (self.OR * chi))
        self.adj = self.Dj.T @ (chi ** 2)


# ---------------------------------------------------------------- helix
def plus(theta):
    return np.array([1.0, np.exp(1j * theta)]) / np.sqrt(2)


def helix_track(N, phi, pin):
    """Full 2^N helix track: sum of h_phi bonds plus pin a=I-|+0><+0| on site `pin`."""
    s = np.zeros(4, complex)
    s[1] = 1 / np.sqrt(2)              # |01>
    s[2] = -np.exp(-1j * phi) / np.sqrt(2)   # |10>
    h = np.outer(s, s.conj())
    a = np.eye(2) - np.outer(plus(0), plus(0).conj())
    dim = 2 ** N
    H = sp.csr_matrix((dim, dim), dtype=complex)
    for x in range(N - 1):
        H = H + sp.kron(sp.kron(sp.identity(2 ** x), h), sp.identity(2 ** (N - x - 2)))
    H = H + sp.kron(sp.kron(sp.identity(2 ** pin), a), sp.identity(2 ** (N - pin - 1)))
    return H.toarray()


def site_op(N, x, o):
    return np.kron(np.kron(np.eye(2 ** x), o), np.eye(2 ** (N - x - 1)))


def compensator(b, phi, r1, r2, kind='small'):
    u1, u2 = plus(-(b - 1) * phi), plus(-b * phi)
    m1, m2 = 1 / r1, 1 / r2
    chi = np.sin(phi / 2) ** 2
    rho1, rho2 = np.outer(u1, u1.conj()), np.outer(u2, u2.conj())
    if kind == 'small':          # C8 (8a)
        M = (m1 + m2) / 2 * np.eye(2) + (m1 - m2) / (2 * chi) * (rho1 - rho2)
    else:                        # C8 (4), C44 (8)
        M = (m1 * (np.eye(2) - rho2) + m2 * (np.eye(2) - rho1)) / chi
    return M, (u1, u2)


# ---------------------------------------------------------------- selector head
def head(l):
    """Sine-path ground state of the selector head on l-1 positions."""
    i = np.arange(1, l)
    return np.sqrt(2 / l) * np.sin(np.pi * i / l)


# ---------------------------------------------------------------- wrapped cell
class WrappedCell:
    """Charge x core(bus r excited) operator H_(p,r) of C44 (PROOFS Sec. 5)."""

    def __init__(self, core, ch, s, bus):
        self.core, self.ch, self.s, self.bus = core, ch, s, bus
        idx = core.sector(bus)
        self.idx = idx
        dc = len(idx)
        bonds = [b.toarray()[np.ix_(idx, idx)] for b in core.bonds]
        ends = (core.aL + core.bR).toarray()[np.ix_(idx, idx)]
        m = len(ch.configs)
        H = np.kron(ch.H - ch.e * np.eye(m), np.eye(dc))
        W = np.kron(np.eye(m), ends)
        for j, bj in enumerate(bonds):
            W = W + np.kron(np.diag(ch.Dj[:, j]), bj)
        self.W = W
        self.Hch = np.kron(ch.H - ch.e * np.eye(m), np.eye(dc))
        self.H = self.Hch + s * W
        self.dim = H.shape[0]

    def source_vectors(self, end, G):
        """dict sigma -> matrix (dim x 2^n) of (O_e chi) x F^e_sigma G|i>."""
        O = self.ch.OL if end == 'A' else self.ch.OR
        out = {}
        for sg in SPECIES:
            F = self.core.source(end, self.bus, sg).toarray()
            cvec = (F @ G)[self.idx]
            out[sg] = np.kron((O * self.ch.chi)[:, None], cvec)
        return out

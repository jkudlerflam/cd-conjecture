"""Open-chain gap representation of the legally projected vacuum-clock vacancy model (C89, C94).

Material sites j = 0..N-1 in order. Site j carries g_j, the number of vacancies between material j
and j+1 (j <= N-2; the last site carries no gap, so vacancies never leave the material block), and
a clock bit c_j. Cells are consecutive runs of material sites with lengths `lengths`; a non-pinned
cell carries exactly one clock, a pinned cell none. Hamiltonian (t = 1 by default):

  H = -t sum_{j=0}^{N-3} (S^-_j S^+_{j+1} + h.c.)                vacancy hopping, norm-one shifts
      -J sum_{cells, q} A_{start+q} (clock hop start+q <-> start+q+1)
      + J mu/2 sum_{markers m} (A_{m-1} + A_m)                     (A_{N-1} absent on the open chain)

with A_j = [g_j = 0] the contact projector. This is C89/C94's `clock_model` with the ring replaced
by an open block of material, all signs and coefficients unchanged.
"""
import math
import itertools
import numpy as np


def selector(b):
    def a(l): return 2 * (1 - math.cos(math.pi / (l + 1)))
    def d(l): return (l + 1) * a(l) - l * a(l + 1)
    return 2 - (d(b) + d(b - 1)) / 2


def cell_energy(l, mu):
    return mu - 2 * math.cos(math.pi / (l + 1))


class Model:
    def __init__(self, lengths, pins, V, J, mu, G, t=1.0):
        self.lengths = list(lengths); self.pins = list(pins)
        self.N = N = sum(lengths); self.V = V; self.J = J; self.mu = mu; self.G = G; self.t = t
        self.d = 2 * (G + 1)
        cell = []; pos = []
        for c, b in enumerate(lengths):
            cell += [c] * b; pos += list(range(b))
        self.cell = cell; self.pos = pos
        self.markers = list(np.cumsum(lengths) - 1)
        self.qloc = np.array([[g, c] for g in range(G + 1) for c in (0, 1)])
        # allowed clock counts at each bond k (cumulative over sites < k)
        self.done = []; self.inside = []
        for k in range(N + 1):
            nd = sum(1 for c in range(len(lengths)) if not pins[c] and sum(lengths[:c + 1]) <= k)
            ins = 0
            if 0 < k < N and cell[k - 1] == cell[k] and not pins[cell[k]]:
                ins = 1
            self.done.append(nd); self.inside.append(ins)

    def s(self, g, c):
        return 2 * g + c

    def allowed(self, k, q):
        v, c = q
        if v < 0 or v > self.V:
            return False
        if k >= self.N - 1 and v != self.V:
            return False
        return self.done[k] <= c <= self.done[k] + self.inside[k]

    def site_diag(self, j):
        """Diagonal single-site energy on site j (coefficient of A_j)."""
        coef = 0.0
        if j <= self.N - 2:
            for m in self.markers:
                if j == m - 1 or j == m:
                    coef += self.J * self.mu / 2
        return coef

    def bonds(self):
        d = self.d; G = self.G; N = self.N
        hs = []
        for j in range(N - 1):
            h = np.zeros((d * d, d * d))
            def ix(s1, s2): return s1 * d + s2
            # vacancy hopping between gaps j and j+1 (both internal)
            if j + 1 <= N - 2:
                for g1 in range(G + 1):
                    for g2 in range(G + 1):
                        for c1 in (0, 1):
                            for c2 in (0, 1):
                                if g1 >= 1 and g2 + 1 <= G:
                                    a = ix(self.s(g1, c1), self.s(g2, c2)); b = ix(self.s(g1 - 1, c1), self.s(g2 + 1, c2))
                                    h[b, a] -= self.t; h[a, b] -= self.t
            # clock hop j <-> j+1 inside a non-pinned cell, weighted by A_j
            if self.cell[j] == self.cell[j + 1] and not self.pins[self.cell[j]]:
                for g2 in range(G + 1):
                    a = ix(self.s(0, 1), self.s(g2, 0)); b = ix(self.s(0, 0), self.s(g2, 1))
                    h[b, a] -= self.J; h[a, b] -= self.J
            # single-site contact energies: site j on bond j; the last site has none
            cj = self.site_diag(j)
            for g2 in range(G + 1):
                for c1 in (0, 1):
                    for c2 in (0, 1):
                        a = ix(self.s(0, c1), self.s(g2, c2))
                        h[a, a] += cj
            hs.append(h)
        return hs

    def config(self):
        """Product state: vacancies spread over the N-1 internal gaps, clocks at cell centres."""
        N, V = self.N, self.V
        gaps = np.zeros(N, int)
        for k in range(V):
            gaps[int((k + 0.5) * (N - 1) / V)] += 1
        assert gaps.max() <= self.G and gaps[N - 1] == 0
        cfg = []
        start = 0
        clock = np.zeros(N, int)
        for c, b in enumerate(self.lengths):
            if not self.pins[c]:
                clock[start + (b - 1) // 2] = 1
            start += b
        return [self.s(int(g), int(c)) for g, c in zip(gaps, clock)]

    def first_order(self, contact):
        """J times the ground energy of the first-order clock Hamiltonian with contact profile
        contact[j] = <A_j>_0 (free bath ground state), j = 0..N-2."""
        A = np.append(np.asarray(contact, float), 0.0)
        E = 0.0
        start = 0
        for c, b in enumerate(self.lengths):
            m = start + b - 1
            E += self.mu / 2 * (A[m - 1] + A[m])
            if not self.pins[c]:
                T = np.zeros((b, b))
                for q in range(b - 1):
                    T[q, q + 1] = T[q + 1, q] = -A[start + q]
                E += np.linalg.eigvalsh(T)[0]
            start += b
        return self.J * E

    def frozen(self):
        return sum(self.mu if p else cell_energy(b, self.mu) for b, p in zip(self.lengths, self.pins))


def free_energy_exact(N, V, t=1.0):
    """J = 0 ground energy: N-2 movable material and V vacancies on the L-2 interior sites with the
    end material fixed, i.e. V free fermions (the vacancies) on an open chain of L-2 sites."""
    L = N + V
    n = L - 2
    k = np.arange(1, n + 1)
    e = np.sort(-2 * t * np.cos(np.pi * k / (n + 1)))
    return float(e[:V].sum())


def ed_open(lengths, pins, V, J, mu, t=1.0):
    """Independent sparse ED in the C89 construction (gap compositions x clock product), open chain."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as sla
    N = sum(lengths); ng = N - 1
    gaps = []
    for sep in itertools.combinations(range(V + ng - 1), ng - 1):
        p = (-1,) + sep + (V + ng - 1,)
        gaps.append(tuple(p[j + 1] - p[j] - 1 for j in range(ng)))
    look = {g: i for i, g in enumerate(gaps)}
    rows, cols = [], []
    for i, g in enumerate(gaps):
        for j in range(ng):
            if g[j]:
                for dest in (j - 1, j + 1):
                    if 0 <= dest < ng:
                        h = list(g); h[j] -= 1; h[dest] += 1
                        rows.append(look[tuple(h)]); cols.append(i)
    dg = len(gaps)
    Hg = sp.csr_matrix((-t * np.ones(len(rows)), (rows, cols)), shape=(dg, dg))
    Gm = np.asarray(gaps)
    A = np.hstack([(Gm == 0).astype(float), np.zeros((dg, 1))])
    dims = [1 if p else b for b, p in zip(lengths, pins)]
    states = list(itertools.product(*(range(x) for x in dims)))
    idx = {x: i for i, x in enumerate(states)}; dc = len(states)
    markers = np.cumsum(lengths) - 1
    diag = mu / 2 * sum(A[:, m - 1] + A[:, m] for m in markers)
    H = sp.kron(Hg + sp.diags(J * diag), sp.eye(dc), format='csr')
    start = 0
    for c, (b, p) in enumerate(zip(lengths, pins)):
        if not p:
            for q in range(b - 1):
                r = []; cc = []
                for s, i in idx.items():
                    if s[c] == q:
                        z = list(s); z[c] += 1; k = idx[tuple(z)]
                        r += [i, k]; cc += [k, i]
                hop = sp.csr_matrix((-np.ones(len(r)), (r, cc)), shape=(dc, dc))
                H = H + J * sp.kron(sp.diags(A[:, start + q]), hop, format='csr')
        start += b
    e = sla.eigsh(H, k=1, which='SA', tol=1e-13)[0][0]
    e0, v0 = sla.eigsh(Hg, k=1, which='SA', tol=1e-13)
    contact = (v0[:, 0] ** 2) @ A[:, :ng]
    return float(e), float(e0[0]), contact, H.shape[0]

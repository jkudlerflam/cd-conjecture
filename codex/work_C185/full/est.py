"""C155 fast central-charge estimator for a translation-invariant qubit ring (Claude Opus agent, 2026-11-07).

Hamiltonian on a ring of L qubits, H = sum_j T^j h T^-j, with h a 2^r x 2^r matrix acting on sites j..j+r-1
(r = 2 for the nearest-neighbour bonds of B_2; r = 3 is used only for calibration models such as two
decoupled chains or the cluster-Ising chain).

Point convention for B_2 (same as contrib/d2_grid/strata.py): v = (Jx, Jy, Jz, Dx, Dy, Dz, fx, fy, fz),
H = sum_j [s_j^T J s_{j+1} + f . s_j], J = diag(Jx, Jy, Jz) + A(D), A_ab = eps_abc D_c. v is normalized
to the unit 8-sphere (principal frame of the symmetric part, which is the SO(3)-reduced parametrization).

Momentum-resolved Lanczos: the sector matrices are assembled from a precomputed table, so that a new bond
only needs a gather and a bincount (no sparse-matrix products).
Basis convention: site i is bit (L-1-i) of the integer (site 0 most significant), as in work_C143/ed.py.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [SX, SY, SZ]


def J_f(v):
    v = np.asarray(v, float)
    J = np.diag(v[:3]).astype(float)
    Dx, Dy, Dz = v[3:6]
    J[0, 1] += Dz; J[1, 0] -= Dz
    J[1, 2] += Dx; J[2, 1] -= Dx
    J[2, 0] += Dy; J[0, 2] -= Dy
    return J, v[6:9].copy()


def bond_from_JB(J, f):
    """4x4 bond, field split symmetrically over the two sites (ring-equivalent to f . s_j per site)."""
    h = np.zeros((4, 4), complex)
    for a in range(3):
        for b in range(3):
            h += J[a, b] * np.kron(PAULI[a], PAULI[b])
        h += f[a] * (np.kron(PAULI[a], I2) + np.kron(I2, PAULI[a])) / 2
    return h


def bond_v(v):
    v = np.asarray(v, float)
    v = v / np.linalg.norm(v)
    return bond_from_JB(*J_f(v))


# ---------------------------------------------------------------- momentum tables
_orb = {}


def orbits(L):
    if L in _orb:
        return _orb[L]
    N = 1 << L
    x = np.arange(N, dtype=np.int64)
    mask = N - 1
    rep = x.copy(); jmin = np.zeros(N, dtype=np.int64)
    cur = x.copy()
    per = np.full(N, L, dtype=np.int64)
    for j in range(1, L):
        cur = ((cur >> 1) | ((cur & 1) << (L - 1))) & mask   # T: site i -> i+1
        better = cur < rep
        rep[better] = cur[better]; jmin[better] = j
        hit = (cur == x) & (per == L)
        per[hit] = j
    # rep = T^jmin x  ->  x = T^{-jmin} rep, shift = -jmin mod L  (x = T^shift rep)
    shift = (-jmin) % L
    _orb[L] = (rep, shift, per)
    return _orb[L]


class RingTable:
    """Precomputed sector tables for ring length L and bond range r."""

    def __init__(self, L, r=2, keep=True):
        self.L, self.r, self.keep = L, r, keep
        rep, shift, per = orbits(L)
        N = 1 << L
        self.N = N
        x = np.arange(N, dtype=np.int64)
        self.reps_all = x[rep == x]
        self.rep, self.shift, self.per = rep, shift, per
        self._sec = {}

    def sector(self, m):
        if m in self._sec:
            return self._sec[m]
        L, r = self.L, self.r
        rep, shift, per, N = self.rep, self.shift, self.per, self.N
        R_all = self.reps_all
        k = 2 * np.pi * m / L
        ok = ((m * per[R_all]) % L) == 0
        R = R_all[ok]
        n = len(R)
        pos = -np.ones(N, dtype=np.int64)
        pos[R] = np.arange(n)
        rows, cols, hidx, coef = [], [], [], []
        dout = 1 << r
        for i in range(L):
            sites = [(i + s) % L for s in range(r)]
            bits = [L - 1 - s for s in sites]
            inidx = np.zeros(n, dtype=np.int64)
            base = R.copy()
            for s, b in enumerate(bits):
                inidx |= ((R >> b) & 1) << (r - 1 - s)
                base &= ~(1 << b)
            for out in range(dout):
                y = base.copy()
                for s, b in enumerate(bits):
                    y |= ((out >> (r - 1 - s)) & 1) << b
                ry = rep[y]
                col_ok = pos[ry] >= 0
                if not col_ok.any():
                    continue
                yy = y[col_ok]
                rows.append(pos[ry[col_ok]])          # <r_y, k| H |r, k>
                cols.append(np.nonzero(col_ok)[0])
                hidx.append(out * dout + inidx[col_ok])
                coef.append(np.exp(1j * k * shift[yy]) * np.sqrt(per[R[col_ok]] / per[ry[col_ok]]))
        rows = np.concatenate(rows); cols = np.concatenate(cols)
        hidx = np.concatenate(hidx).astype(np.int16); coef = np.concatenate(coef)
        key = rows * n + cols
        ukey, inv = np.unique(key, return_inverse=True)
        urow = ukey // n; ucol = ukey % n
        indptr = np.zeros(n + 1, dtype=np.int64)
        np.add.at(indptr, urow + 1, 1)
        indptr = np.cumsum(indptr)
        sec = dict(m=m, k=k, R=R, n=n, inv=inv.astype(np.int32), hidx=hidx, coef=coef,
                   indices=ucol.astype(np.int32), indptr=indptr, nnz=len(ukey))
        if self.keep:
            self._sec[m] = sec
        return sec

    def sector_matrix(self, h, m):
        s = self.sector(m)
        hv = h.reshape(-1)[s['hidx']] * s['coef']
        data = np.bincount(s['inv'], weights=hv.real, minlength=s['nnz']) + \
            1j * np.bincount(s['inv'], weights=hv.imag, minlength=s['nnz'])
        return sp.csr_matrix((data, s['indices'], s['indptr']), shape=(s['n'], s['n']))

    def full_vector(self, phi, m):
        """Momentum-basis vector -> computational basis."""
        s = self.sector(m)
        rep, shift, per = self.rep, self.shift, self.per
        pos = -np.ones(self.N, dtype=np.int64)
        pos[s['R']] = np.arange(s['n'])
        p = pos[rep]
        ok = p >= 0
        psi = np.zeros(self.N, complex)
        xs = np.nonzero(ok)[0]
        psi[xs] = phi[p[xs]] * np.exp(-1j * s['k'] * shift[xs]) / np.sqrt(per[xs])
        return psi / np.linalg.norm(psi)


_tables = {}


def table(L, r=2):
    if (L, r) not in _tables:
        _tables[(L, r)] = RingTable(L, r, keep=(L <= 16))
    return _tables[(L, r)]


def lowest(Hk, nev):
    n = Hk.shape[0]
    if n <= 400:
        w, V = np.linalg.eigh(Hk.toarray())
        return w[:nev], V[:, :nev]
    w, V = spla.eigsh(Hk, k=nev, which='SA', tol=1e-10, ncv=max(20, 2 * nev + 8))
    o = np.argsort(w)
    return w[o], V[:, o]


def ring_data(h, L, r=2, nev=2):
    """Lowest nev levels in every momentum sector, and the ground state (full basis)."""
    T = table(L, r)
    lev = np.full((L, nev), np.inf)
    best = (np.inf, None, None)
    for m in range(L):
        Hk = T.sector_matrix(h, m)
        if Hk.shape[0] == 0:
            continue
        ne = min(nev, Hk.shape[0])
        w, V = lowest(Hk, ne)
        lev[m, :ne] = w
        if w[0] < best[0] - 1e-12 or (abs(w[0] - best[0]) <= 1e-12 and min(m, L - m) < min(best[1], L - best[1])):
            best = (w[0], m, V[:, 0])
    E0, m0, phi = best
    psi = T.full_vector(phi, m0)
    return dict(L=L, lev=lev, E0=E0, m0=m0, psi=psi)


def entropies(psi, L):
    """Von Neumann entropy of sites 0..l-1 for l = 1..L/2."""
    S = np.zeros(L // 2 + 1)
    for l in range(1, L // 2 + 1):
        s = np.linalg.svd(psi.reshape(1 << l, 1 << (L - l)), compute_uv=False)
        p = s ** 2; p = p[p > 1e-14]; p /= p.sum()
        S[l] = -np.sum(p * np.log(p))
    return S


def estimate(h, Ls=(12, 16), r=2, nev=2):
    """Returns the screening estimators.

    c_S   : 3 [S_{L2}(L2/2) - S_{L1}(L1/2)] / log(L2/L1), half-ring entropies (velocity free).
    c_fit : slope of S_L(l) against (1/3) log[(L/pi) sin(pi l/L)] for l = 2..L/2 at the largest L.
    c_cas : Casimir product (pi v c / 6) from E0/L at the two sizes, divided by v from the lowest level at
            momentum m0 +- 1 at the largest size (v = L dE/2pi, scaling dimension 1 assumed there).
    c_obj : min(c_S, c_cas), the quantity maximized by the search (both must be fooled together); set to 0
            when three or more levels lie within 1e-3/L of the ground state at L2 (degenerate manifold,
            e.g. the XZ point or FM multiplets, where the ground-state entropy is not defined).
    """
    out = {}
    datas = [ring_data(h, L, r, nev) for L in Ls]
    L1, L2 = Ls[0], Ls[-1]
    d1, d2 = datas[0], datas[-1]
    S1 = entropies(d1['psi'], L1); S2 = entropies(d2['psi'], L2)
    c_S = 3 * (S2[L2 // 2] - S1[L1 // 2]) / np.log(L2 / L1)
    ls = np.arange(2, L2 // 2 + 1)
    X = np.log(L2 / np.pi * np.sin(np.pi * ls / L2)) / 3
    c_fit = np.polyfit(X, S2[ls], 1)[0]
    e1, e2 = d1['E0'] / L1, d2['E0'] / L2
    vc = (6 / np.pi) * (e2 - e1) / (1 / L1 ** 2 - 1 / L2 ** 2)
    m0 = d2['m0']
    lev = d2['lev']
    cand = [lev[(m0 + 1) % L2, 0], lev[(m0 - 1) % L2, 0]]
    dE = min(cand) - d2['E0']
    v = L2 * dE / (2 * np.pi)
    # a vanishing velocity means a degenerate level at m0 +- 1, where the Casimir ratio is meaningless
    c_cas = vc / v if v > 0.02 else 0.0
    # gap above the ground state (any momentum) at both sizes, scaled by L
    gaps = []; ndeg = []; gap2 = []
    for d in datas:
        w = np.sort(d['lev'].ravel())
        gaps.append(float(w[1] - w[0])); gap2.append(float(w[2] - w[0]))
        ndeg.append(int(np.sum(d['L'] * (w - w[0]) < 1e-3)))
    # effective dynamical exponents from the first and second gaps, z = -dlog gap / dlog L
    zs = []
    for g in (gaps, gap2):
        zs.append(float(-np.log(max(g[-1], 1e-14) / max(g[0], 1e-14)) / np.log(L2 / L1)))
    out['z1'], out['z2'] = zs
    # criticality indicator: |z1 - 1|, or |z2 - 1| when the first gap closes faster than 1/L^1.8
    # (exponentially split cat or dimer partner, so the second gap is the physical one)
    out['crit'] = float(abs(zs[0] - 1) if zs[0] <= 1.8 else abs(zs[1] - 1))
    out.update(c_S=float(c_S), c_fit=float(c_fit), c_cas=float(c_cas), vc=float(vc), v=float(v),
               c_obj=float(min(c_S, c_cas)) if ndeg[-1] < 3 else 0.0, ndeg=ndeg, gap1=gaps, gap2=gap2, E0=[float(d['E0']) for d in datas],
               m0=[int(d['m0']) for d in datas], Lgap=[float(g * d['L']) for g, d in zip(gaps, datas)],
               Shalf=[float(S1[L1 // 2]), float(S2[L2 // 2])])
    return out


# ----------------------------------------------------------------- self tests
if __name__ == '__main__':
    import time, sys
    sys.path.insert(0, '../work_C143')
    import ed as ed143
    rng = np.random.default_rng(1)
    for L in (6, 8):
        v = rng.normal(size=9)
        h = bond_v(v)
        T = table(L)
        H = ed143.ham_sparse(h, L).toarray()
        wd = np.linalg.eigvalsh(H)
        allw = []
        for m in range(L):
            Hk = T.sector_matrix(h, m).toarray()
            assert np.allclose(Hk, Hk.conj().T), 'non-Hermitian sector'
            allw.append(np.linalg.eigvalsh(Hk))
            # compare with ed143 projector
            P = ed143.momentum_projector(L, m)
            Hk2 = (P.conj().T @ ed143.ham_sparse(h, L) @ P).toarray()
            assert np.allclose(np.sort(np.linalg.eigvalsh(Hk)), np.sort(np.linalg.eigvalsh(Hk2)))
        allw = np.sort(np.concatenate(allw))
        print('L', L, 'spectrum diff vs dense', np.max(np.abs(allw - wd)))
        dd = ring_data(h, L)
        psi = dd['psi']
        print('  ground-state residual', np.linalg.norm(H @ psi - dd['E0'] * psi), 'E0 diff', dd['E0'] - wd[0])
    # range-3 check
    L = 8
    h3 = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8)); h3 = h3 + h3.conj().T
    T = table(L, 3)
    w = np.sort(np.concatenate([np.linalg.eigvalsh(T.sector_matrix(h3, m).toarray()) for m in range(L)]))
    # dense range-3 ring
    H = np.zeros((1 << L, 1 << L), complex)
    for i in range(L):
        sites = [(i + s) % L for s in range(3)]
        # permute: build operator by acting on basis states
        for x in range(1 << L):
            bits = [(x >> (L - 1 - s)) & 1 for s in sites]
            inn = bits[0] * 4 + bits[1] * 2 + bits[2]
            base = x
            for s in sites:
                base &= ~(1 << (L - 1 - s))
            for out in range(8):
                y = base
                for t, s in enumerate(sites):
                    y |= ((out >> (2 - t)) & 1) << (L - 1 - s)
                H[y, x] += h3[out, inn]
    print('range-3 L=8 diff', np.max(np.abs(w - np.linalg.eigvalsh(H))))
    for Ls in [(12, 16)]:
        t = time.time()
        for L in Ls:
            table(L)
        print('tables', Ls, time.time() - t)
        t = time.time()
        print(estimate(bond_v(rng.normal(size=9)), Ls))
        print('one estimate', time.time() - t)

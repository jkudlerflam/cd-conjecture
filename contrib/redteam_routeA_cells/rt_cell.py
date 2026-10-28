"""Cells and chains generated from the literal 23-letter bond of rt_bond.py.

Legal words are enumerated from the arrow grammar (not from a hand-written
cell list), every Hamiltonian is assembled by applying the literal bond on
each nearest-neighbour pair, and closure of the enumerated space is asserted.
The second-order bond is obtained by a sparse direct solve in the two-cell
excited space (not by the Cauchy eigen-sum used in contrib/routeA_spatial/).
"""
import itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from rt_bond import LETTERS, IX, D0, arrows, legal_pair, charge, is_exc, is_marker, bond, onsite_matrix, F_ops, P

NONM = [i for i, x in enumerate(LETTERS) if not is_marker(x)]
MARK = [i for i, x in enumerate(LETTERS) if is_marker(x)]
BMARK = IX[('M', None, None)]


def bond_dict(h):
    """h as {(i,j): [(coef,(i2,j2)), ...]} for the off-diagonal and diagonal action."""
    out = {}
    nz = np.nonzero(np.abs(h) > 0)
    for r, c in zip(*nz):
        out.setdefault((c // D0, c % D0), []).append((h[r, c], (r // D0, r % D0)))
    return out


def enum_words(nsites, marker_pos, left_fixed=BMARK, cyclic=False):
    """All words on nsites sites with markers exactly at marker_pos and no
    forbidden pair (including the pair with a fixed left boundary letter, or the
    wrap-around pair if cyclic)."""
    words = []

    def rec(w):
        k = len(w)
        if k == nsites:
            if cyclic and not legal_pair(LETTERS[w[-1]], LETTERS[w[0]]):
                return
            words.append(tuple(w))
            return
        cand = MARK if k in marker_pos else NONM
        prev = w[-1] if k > 0 else left_fixed
        for c in cand:
            if prev is not None and not legal_pair(LETTERS[prev], LETTERS[c]):
                continue
            rec(w + [c])
    rec([])
    return words


def build(words, bonds, hdict, onsite_sites, h1, left_fixed=None, check_closure=True):
    """Sparse H on span(words).  bonds: list of (i, j, hdict_for_that_bond); site -1
    denotes the fixed left boundary letter.  onsite_sites: sites receiving h1."""
    ix = {w: k for k, w in enumerate(words)}
    rows, cols, vals = [], [], []
    leak = 0.0
    h1d = {}
    for i in range(D0):
        for j in range(D0):
            if h1[j, i] != 0:
                h1d.setdefault(i, []).append((h1[j, i], j))
    for k, w in enumerate(words):
        for (s1, s2, hd) in bonds:
            a = left_fixed if s1 == -1 else w[s1]
            b = w[s2]
            for c, (a2, b2) in hd.get((a, b), []):
                if s1 == -1 and a2 != left_fixed:
                    leak += abs(c) ** 2
                    continue
                w2 = list(w)
                if s1 != -1:
                    w2[s1] = a2
                w2[s2] = b2
                w2 = tuple(w2)
                if w2 not in ix:
                    leak += abs(c) ** 2
                    continue
                rows.append(ix[w2]); cols.append(k); vals.append(c)
        for s in onsite_sites:
            for c, b in h1d.get(w[s], []):
                w2 = list(w); w2[s] = b; w2 = tuple(w2)
                rows.append(ix[w2]); cols.append(k); vals.append(c)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(len(words), len(words)))
    if check_closure:
        assert leak < 1e-20, ('basis not closed', leak)
    return H


def cell(n, p):
    """One cell: n data sites + its marker, left neighbour = frozen marker letter
    (supplies the junction).  Bonds carry kinetic+penalty only (no port inside a
    cell); full on-site terms on the n+1 cell sites."""
    hk = bond(p, lam=0.0, onsite=False)
    hd = bond_dict(hk)
    words = enum_words(n + 1, {n}, left_fixed=BMARK)
    bonds = [(-1, 0, hd)] + [(k, k + 1, hd) for k in range(n)]
    H = build(words, bonds, hd, range(n + 1), onsite_matrix(p), left_fixed=BMARK)
    Hd = H.toarray()
    assert np.abs(Hd - Hd.conj().T).max() < 1e-13
    info = dict(words=words, Q=np.array([sum(charge(LETTERS[c]) for c in w) for w in words]),
                exc=np.array([any(is_exc(LETTERS[c]) for c in w) for w in words]),
                data=[tuple(LETTERS[c][2] for c in w[:n]) for w in words])
    return Hd, info


def code_basis(n, Hd, info, tol=1e-10):
    """Ground space, its degeneracy, gap, and the canonical basis |walker> x |a>
    (data configuration a), built by projecting a reference word."""
    w, U = np.linalg.eigh(Hd)
    N = 2 ** n
    E0 = w[0]
    deg = int(np.sum(w < E0 + tol))
    gap = w[deg] - E0
    Pg = U[:, :deg] @ U[:, :deg].conj().T
    ix = {wd: k for k, wd in enumerate(info['words'])}
    G = np.zeros((len(info['words']), N), complex)
    for col in range(N):
        a = tuple(int(b) for b in format(col, '0%db' % n)) if n else ()
        # reference: vacuum token on the first data site, data a
        ref = (IX[('D', 'v', a[0])],) + tuple(IX[('r', None, a[j])] for j in range(1, n)) + (IX[('M', None, None)],)
        v = Pg[:, ix[ref]]
        v = v / np.linalg.norm(v)
        v = v * np.exp(-1j * np.angle(v[ix[ref]]))
        G[:, col] = v
    orth = float(np.abs(G.conj().T @ G - np.eye(N)).max())
    # each code vector has definite physical data and vacuum token
    purity = 0.0
    for col in range(N):
        a = tuple(int(b) for b in format(col, '0%db' % n)) if n else ()
        wrong = sum(abs(G[k, col]) ** 2 for k in range(len(info['words']))
                    if info['data'][k] != a or info['exc'][k])
        purity = max(purity, wrong)
    return dict(E0=E0, deg=deg, gap=gap, G=G, orth_err=orth, impurity=purity,
                eig_spread=float(w[deg - 1] - w[0]), evals=w)


def port_factors(n, info):
    """F^M on the cell's own marker (right port), F^X on its first site (left port)."""
    F = F_ops()
    words = info['words']
    ix = {w: k for k, w in enumerate(words)}
    FM, FX = {}, {}
    for f in F:
        A = np.zeros((len(words), len(words)))
        B = np.zeros((len(words), len(words)))
        for k, w in enumerate(words):
            for (pos, store, req) in ((n, A, 'M'), (0, B, 'D')):
                c = w[pos]
                if LETTERS[c][0] != req:
                    continue
                col = F[f][:, c]
                for c2 in np.nonzero(col)[0]:
                    w2 = list(w); w2[pos] = c2; w2 = tuple(w2)
                    store[ix[w2], k] += col[c2]
        FM[f], FX[f] = A, B
    return FM, FX


def port_V(p, FM, FX):
    """lam=1 inter-cell bond on (left cell marker) x (right cell first site)."""
    return [(FM['0'], FX['0'], 1.0), (FM['p'], FX['m'], np.exp(1j * p['phi'])),
            (FM['m'], FX['p'], np.exp(-1j * p['phi']))]


def second_order(n, p, return_all=False):
    Hd, info = cell(n, p)
    cb = code_basis(n, Hd, info)
    G, E0 = cb['G'], cb['E0']
    FM, FX = port_factors(n, info)
    terms = port_V(p, FM, FX)
    ex = info['exc']
    ie = np.where(ex)[0]
    He = sp.csr_matrix(Hd[np.ix_(ie, ie)])
    De = len(ie)
    Ie = sp.identity(De, format='csr')
    K = (sp.kron(He, Ie) + sp.kron(Ie, He) - 2 * E0 * sp.identity(De * De)).tocsc()
    N = G.shape[1]
    # X = V (G x G), restricted to exc x exc; check nothing lands elsewhere
    X = np.zeros((De * De, N * N), complex)
    outside = 0.0
    for (A, B, c) in terms:
        AG, BG = A @ G, B @ G
        outside += np.linalg.norm(AG[~ex]) ** 2 + np.linalg.norm(BG[~ex]) ** 2
        X += c * np.einsum('ia,jb->ijab', AG[ie], BG[ie]).reshape(De * De, N * N)
    lu = spla.splu(K)
    Y = lu.solve(X)
    H2 = -(X.conj().T @ Y)
    # first-order compressions G^dag F G
    first = max(float(np.abs(G.conj().T @ M @ G).max()) for M in list(FM.values()) + list(FX.values()))
    res = dict(n=n, E0=E0, deg=cb['deg'], gap=cb['gap'], code_orth=cb['orth_err'],
               code_impurity=cb['impurity'], code_spread=cb['eig_spread'], first_order=first,
               V_leaves_exc=outside, H2=H2, dim_cell=len(info['words']))
    if return_all:
        res.update(Hd=Hd, info=info, G=G, FM=FM, FX=FX, terms=terms)
    return res


# ------------------------------------------------------------- fermion tools
def jw(m):
    Z = np.diag([1.0, -1.0]); I2 = np.eye(2); sm = np.array([[0, 1], [0, 0.]])  # bit 1 occupied
    out = []
    for j in range(m):
        M = np.array([[1.0]])
        for k in range(m):
            M = np.kron(M, Z if k < j else (sm if k == j else I2))
        out.append(M)
    return out


def one_particle_T(H2, n):
    m = 2 * n
    idx = [1 << (m - 1 - j) for j in range(m)]   # mode j <-> bit j in big-endian order
    T = np.array([[H2[idx[x], idx[n + y]] for y in range(n)] for x in range(n)])
    muA = np.array([[H2[idx[x], idx[y]] for y in range(n)] for x in range(n)]) - H2[0, 0] * np.eye(n)
    muB = np.array([[H2[idx[n + x], idx[n + y]] for y in range(n)] for x in range(n)]) - H2[0, 0] * np.eye(n)
    return T, muA, muB


def hs_quadratic(H2, n):
    """Hilbert-Schmidt projection on span{1, f_i^dag f_j}; returns e0, h, residual norm."""
    m = 2 * n
    c = jw(m)
    ops = [np.eye(2 ** m)] + [c[i].T @ c[j] for i in range(m) for j in range(m)]
    A = np.array([o.reshape(-1) for o in ops]).T
    coef, *_ = np.linalg.lstsq(A, H2.reshape(-1), rcond=None)
    R = H2 - (A @ coef).reshape(H2.shape)
    Ntot = sum(ci.T @ ci for ci in c)
    return coef[0], coef[1:].reshape(m, m), float(np.linalg.norm(R, 2)), float(np.abs(H2 @ Ntot - Ntot @ H2).max())


def herm_part_analysis(T):
    """Eigenvalues of the Hermitian part of e^{-i Phi}T at Phi = arg T00, and the
    best rotation (Crawford number of T / |T|)."""
    sc = np.abs(T).max()
    R = T * np.exp(-1j * np.angle(T[0, 0]))
    e_fixed = np.linalg.eigvalsh((R + R.conj().T) / 2) / sc
    best = -np.inf
    for th in np.linspace(0, 2 * np.pi, 3601):
        Rt = T * np.exp(-1j * th)
        best = max(best, np.linalg.eigvalsh((Rt + Rt.conj().T) / 2)[0] / sc)
    asym = float(np.linalg.norm(R - R.T) / np.linalg.norm(R))
    return e_fixed, best, asym


def bands(T, mu, nk=200001):
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    n = T.shape[0]
    E = np.empty((nk, n))
    for i in range(0, nk, 5000):
        kk = ks[i:i + 5000]
        Hk = mu[None] + T[None] * np.exp(1j * kk)[:, None, None] + T.conj().T[None] * np.exp(-1j * kk)[:, None, None]
        E[i:i + 5000] = np.linalg.eigvalsh(Hk)
    out = []
    for b in range(n):
        e = E[:, b]
        s = np.sign(e)
        ch = np.where(s != np.roll(s, -1))[0]
        roots = [ks[i] - e[i] * (2 * np.pi / nk) / (e[(i + 1) % nk] - e[i]) for i in ch]
        occ = ks[e < 0]
        centre = float(np.angle(np.exp(1j * occ).sum())) if 0 < len(occ) < nk else None
        out.append(dict(roots=roots, centre=centre, filling=len(occ) / nk))
    return out

"""Folded-bus (Y) spatial cell for route A.  numpy/scipy only.

Alphabet (23 letters): L_a, R_a (4); Xd_{s,a}, s in v,0,+,- (8);
Xx_{s,a}, s in 0,+,- (6); M_empty, M_s (5).  Xd is the token on the data
track (stub), Xx the token on the express track.  The vacuum species v only
lives on the data track and the marker, so the selector path is C22's.
Excited species reach the marker only through the express track, and switch
between tracks only on the site to the right of a marker (the junction).
"""
import itertools, json, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

SPECIES = ['v', '0', '+', '-']
EXC = ['0', '+', '-']
CHG = {'v': 0, '0': 0, '+': 1, '-': -1}

LET = ([('L', None, a) for a in (0, 1)] + [('R', None, a) for a in (0, 1)]
       + [('Xd', s, a) for s in SPECIES for a in (0, 1)]
       + [('Xx', s, a) for s in EXC for a in (0, 1)]
       + [('Me', None, None)] + [('M', s, None) for s in SPECIES])
IDX = {l: i for i, l in enumerate(LET)}
D0 = len(LET)  # 23

DEF = dict(t=1.0, tx=1.0, tsw=1.0, Delta=2.0, Dx=3.0, DM=3.0, delta=0.3,
           gp=0.3, gm=0.3, phs={'0': 0.0, '+': 0.0, '-': 0.0},
           phx={'0': 0.4, '+': 0.0, '-': 0.0}, phsw={'0': 0.0, '+': 0.0, '-': 0.0},
           bphi=np.pi / 2, eps=0.0, musel=0.0,
           td=0.0, phd={'0': 0.0, '+': 0.0, '-': 0.0})


def params(**kw):
    P = {k: (dict(v) if isinstance(v, dict) else v) for k, v in DEF.items()}
    for k, v in kw.items():
        P[k] = dict(v) if isinstance(v, dict) else v
    return P


# ---------------------------------------------------------------- literal bond
def bond_moves(P):
    """Two-site kinetic part of the repeated bond, as {(l,r): [(c,(l',r'))]}.
    Hermitian conjugates are included."""
    B = {}

    def add(l, r, l2, r2, c):
        B.setdefault((IDX[l], IDX[r]), []).append((c, (IDX[l2], IDX[r2])))
        B.setdefault((IDX[l2], IDX[r2]), []).append((np.conj(c), (IDX[l], IDX[r])))
    t = P['t']
    for a in (0, 1):
        for b in (0, 1):
            for s in SPECIES:        # data-track hop, token moves right
                ph = 0.0 if s == 'v' else P['phs'][s]
                z = (-1) ** a if s == '0' else 1     # neutral carries Z of departure
                add(('Xd', s, a), ('R', None, b), ('L', None, a), ('Xd', s, b),
                    -t * np.exp(1j * ph) * z)
            for s in EXC:            # express hop
                z = (-1) ** a if s == '0' else 1
                add(('Xx', s, a), ('R', None, b), ('L', None, a), ('Xx', s, b),
                    -P['tx'] * np.exp(1j * P['phx'][s]) * z)
        for s in EXC:                # express -> marker (excited only)
            z = (-1) ** a if s == '0' else 1
            add(('Xx', s, a), ('Me', None, None), ('L', None, a), ('M', s, None),
                -P['tx'] * np.exp(1j * P['phx'][s]) * z)
        add(('Xd', 'v', a), ('Me', None, None), ('L', None, a), ('M', 'v', None), -t)
        for s in EXC:                # weak direct stub -> marker link (closes the loop)
            if P['td'] != 0:
                z = (-1) ** a if s == '0' else 1
                add(('Xd', s, a), ('Me', None, None), ('L', None, a), ('M', s, None),
                    -P['td'] * np.exp(1j * P['phd'][s]) * z)
    for m in [('Me', None, None)] + [('M', s, None) for s in SPECIES]:
        for s in EXC:                # junction switch on (marker, site 0)
            for a in (0, 1):
                add(m, ('Xd', s, a), m, ('Xx', s, a),
                    -P['tsw'] * np.exp(1j * P['phsw'][s]))
    return B


FMAP = {'0': {'v': '0', '0': 'v'}, '+': {'v': '+', '-': 'v'}, '-': {'+': 'v', 'v': '-'}}


def port_terms(P):
    """V = F0^M F0^X + e^{i phi} F+^M F-^X + e^{-i phi} F-^M F+^X on (marker, Xd)."""
    return [('0', '0', 1.0), ('+', '-', np.exp(1j * P['bphi'])),
            ('-', '+', np.exp(-1j * P['bphi']))]


def onsite(P):
    diag = np.zeros(D0)
    for i, (k, s, a) in enumerate(LET):
        e = 0.0
        if k in ('L', 'R', 'Xd', 'Xx'):
            e += P['eps'] * a
        if s in EXC:
            e += P['Delta'] + {'+': P['delta'], '-': -P['delta'], '0': 0.0}[s]
            e += P['Dx'] if k == 'Xx' else (P['DM'] if k == 'M' else 0.0)
        if k in ('Me', 'M'):
            e += P['musel']
        diag[i] = e
    off = {}

    def add(i, j, c):
        off.setdefault(IDX[i], []).append((c, IDX[j]))
        off.setdefault(IDX[j], []).append((np.conj(c), IDX[i]))
    add(('Xd', '0', 1), ('Xd', '+', 0), P['gp'])   # neutral + data 1 -> plus + data 0
    add(('Xd', '0', 0), ('Xd', '-', 1), P['gm'])   # neutral + data 0 -> minus + data 1
    return diag, off


# ---------------------------------------------------------------- one cell
def cell_states(n):
    """Legal cell words: n data sites then the marker.  Returns letter tuples."""
    out = []
    for pos in [('d', j) for j in range(n)] + [('x', j) for j in range(n)] + [('M', n)]:
        for s in (SPECIES if pos[0] != 'x' else EXC):
            for a in itertools.product((0, 1), repeat=n):
                w = []
                for j in range(n):
                    if pos[0] in 'dx' and j == pos[1]:
                        w.append(IDX[('Xd' if pos[0] == 'd' else 'Xx', s, a[j])])
                    elif j < pos[1]:
                        w.append(IDX[('L', None, a[j])])
                    else:
                        w.append(IDX[('R', None, a[j])])
                w.append(IDX[('M', s, None)] if pos[0] == 'M' else IDX[('Me', None, None)])
                out.append(tuple(w))
    return out


def word_info(w):
    q, exc, data = 0, False, []
    for l in w:
        k, s, a = LET[l]
        if a is not None:
            q += a
            data.append(a)
        if s is not None:
            q += CHG[s]
            exc = exc or (s in EXC)
    return q, exc, tuple(data)


def build_cell(n, P):
    """Cell Hamiltonian from the literal bond, with a frozen empty marker on
    the left supplying the junction (the port term is not used there)."""
    S = cell_states(n)
    ix = {w: i for i, w in enumerate(S)}
    Dm = len(S)
    B = bond_moves(P)
    diag, off = onsite(P)
    H = np.zeros((Dm, Dm), complex)
    left = IDX[('Me', None, None)]
    for i, w in enumerate(S):
        full = (left,) + w
        H[i, i] += sum(diag[l] for l in w)
        for p, l in enumerate(w):
            for c, l2 in off.get(l, []):
                w2 = list(w); w2[p] = l2
                H[ix[tuple(w2)], i] += c
        for p in range(len(full) - 1):
            for c, (l2, r2) in B.get((full[p], full[p + 1]), []):
                f2 = list(full); f2[p] = l2; f2[p + 1] = r2
                assert f2[0] == left
                H[ix[tuple(f2[1:])], i] += c
    assert np.allclose(H, H.conj().T)
    FM, FX = {}, {}
    for f in FMAP:
        A = np.zeros((Dm, Dm), complex); Bx = np.zeros((Dm, Dm), complex)
        for i, w in enumerate(S):
            k, s, a = LET[w[-1]]
            if k == 'M' and s in FMAP[f]:
                w2 = list(w); w2[-1] = IDX[('M', FMAP[f][s], None)]
                A[ix[tuple(w2)], i] += 1
            k, s, a = LET[w[0]]
            if k == 'Xd' and s in FMAP[f]:
                w2 = list(w); w2[0] = IDX[('Xd', FMAP[f][s], a)]
                Bx[ix[tuple(w2)], i] += 1
        FM[f], FX[f] = A, Bx
    info = [word_info(w) for w in S]
    Q = np.array([x[0] for x in info]); exc = np.array([x[1] for x in info])
    return dict(n=n, S=S, H=H, FM=FM, FX=FX, Q=Q, exc=exc,
                data=[x[2] for x in info], terms=port_terms(P), P=P)


def code(cell):
    """Vacuum-token code: returns E0, code isometry G (D x 2^n), gap."""
    n, H = cell['n'], cell['H']
    w, U = np.linalg.eigh(H)
    N = 2 ** n
    E0 = w[0]
    assert np.allclose(w[:N], E0, atol=1e-10), w[:N + 1]
    gap = w[N] - E0
    # canonical basis: data configuration a, vacuum token in its ground state
    vac = ~cell['exc']
    Dm = len(cell['S'])
    G = np.zeros((Dm, N), complex)
    Hv = H[np.ix_(vac, vac)]
    wv, Uv = np.linalg.eigh(Hv)
    idx_v = np.where(vac)[0]
    for col in range(N):
        a = tuple(int(b) for b in format(col, '0%db' % n)) if n else ()
        sel = [k for k, i in enumerate(idx_v) if cell['data'][i] == a]
        vec = np.zeros(Dm, complex)
        sub = Hv[np.ix_(sel, sel)]
        ws, us = np.linalg.eigh(sub)
        v = us[:, 0] * np.sign(us[np.argmax(abs(us[:, 0])), 0])
        vec[idx_v[sel]] = v
        G[:, col] = vec
        assert abs(ws[0] - E0) < 1e-10
    assert np.allclose(H @ G, E0 * G, atol=1e-10)
    return E0, G, gap


# ---------------------------------------------------------------- second order
def second_order(cell, G, E0):
    """Exact lambda^2 coefficient of the bond (left cell's marker, right cell's
    site 0), on code x code, using the Cauchy form over excited eigenstates."""
    H, ex = cell['H'], cell['exc']
    He = H[np.ix_(ex, ex)]
    w, U = np.linalg.eigh(He)
    w = w - E0
    C = 1.0 / (w[:, None] + w[None, :])
    N = G.shape[1]
    A = {f: U.conj().T @ (cell['FM'][f] @ G)[ex] for f in FMAP}
    Bv = {f: U.conj().T @ (cell['FX'][f] @ G)[ex] for f in FMAP}
    H2 = np.zeros((N * N, N * N), complex)
    for (fm, fx, c) in cell['terms']:
        for (fm2, fx2, c2) in cell['terms']:
            a = np.einsum('ei,ej->eij', A[fm].conj(), A[fm2]).reshape(len(w), -1)
            b = np.einsum('ek,el->ekl', Bv[fx].conj(), Bv[fx2]).reshape(len(w), -1)
            Mx = (a.T @ C @ b).reshape(N, N, N, N)       # i j k l
            H2 -= np.conj(c) * c2 * Mx.transpose(0, 2, 1, 3).reshape(N * N, N * N)
    # poles: distinct excited energies with source weight
    wt = sum(np.linalg.norm(A[f], axis=1) ** 2 + np.linalg.norm(Bv[f], axis=1) ** 2 for f in FMAP)
    acc = np.sort(w[wt > 1e-10])
    k = 1 + int(np.sum(np.diff(acc) > 1e-8)) if len(acc) else 0
    return H2, k


# ---------------------------------------------------------------- fermions
def jw_ops(m):
    """Annihilators c_0..c_{m-1} on m qubits, bit 1 = occupied, mode 0 first."""
    Z = np.diag([1.0, -1.0]); I = np.eye(2); sm = np.array([[0, 1], [0, 0.]])
    cs = []
    for j in range(m):
        ops = [Z] * j + [sm] + [I] * (m - j - 1)
        M = ops[0]
        for o in ops[1:]:
            M = np.kron(M, o)
        cs.append(M)
    return cs


def fit_quadratic(Hm, m):
    """Fit Hm = e0 I + sum_ij h_ij c_i^dag c_j + residual.  Returns e0, h, residual."""
    cs = jw_ops(m)
    ops = [np.eye(2 ** m)] + [cs[i].T @ cs[j] for i in range(m) for j in range(m)]
    A = np.array([o.reshape(-1) for o in ops]).T
    coef, *_ = np.linalg.lstsq(A, Hm.reshape(-1), rcond=None)
    fit = (A @ coef).reshape(Hm.shape)
    res = np.linalg.norm(Hm - fit, 2)
    return coef[0], coef[1:].reshape(m, m), res


def logical(H2, n):
    e0, h, res = fit_quadratic(H2, 2 * n)
    muA, muB = h[:n, :n], h[n:, n:]
    T = h[:n, n:]          # coefficient of c^dag_{A,x} c_{B,y}; A is the left cell
    assert np.allclose(h[n:, :n], T.conj().T, atol=1e-12 * max(1, abs(h).max()))
    return dict(e0=e0, muA=muA, muB=muB, T=T, res=res)


def one_particle_T(H2, n):
    """T_xy = <A_x|H2|B_y> in the one-particle sector (mode m <-> index
    2^(2n-1-m)).  Off-diagonal in the cells, so no scalar is subtracted."""
    m = 2 * n
    idx = [2 ** (m - 1 - j) for j in range(m)]
    return np.array([[H2[idx[x], idx[n + y]] for y in range(n)] for x in range(n)])


def bands(T, mu, nk=40001):
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    E = np.array([np.linalg.eigvalsh(mu + T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k))
                  for k in ks])
    out = []
    for b in range(E.shape[1]):
        e = E[:, b]
        sgn = np.sign(e)
        ch = np.where(sgn != np.roll(sgn, -1))[0]
        roots = []
        for i in ch:
            j = (i + 1) % nk
            k0, k1 = ks[i], ks[i] + 2 * np.pi / nk
            roots.append(k0 - e[i] * (k1 - k0) / (e[j] - e[i]))
        occ = ks[e < 0]
        centre = np.angle(np.exp(1j * occ).sum()) if len(occ) and len(occ) < nk else None
        out.append(dict(min=e.min(), max=e.max(), roots=roots, centre=centre))
    return out


def resonance_margin(roots):
    """Min |p1+p2-p3-p4| mod 2pi over Fermi points with {1,2} != {3,4} as multisets,
    and min |p1-p2| (distinct points)."""
    p = np.array(roots)
    def d(x):
        return abs((x + np.pi) % (2 * np.pi) - np.pi)
    m4 = np.inf
    idx = range(len(p))
    pairs = list(itertools.combinations_with_replacement(idx, 2))
    for (a, b) in pairs:
        for (c, e) in pairs:
            if (a, b) != (c, e):
                m4 = min(m4, d(p[a] + p[b] - p[c] - p[e]))
    m2 = min(d(p[a] - p[b]) for a, b in itertools.combinations(idx, 2)) if len(p) > 1 else np.inf
    return m4, m2

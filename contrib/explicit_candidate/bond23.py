"""Independent re-implementation of the 23-letter route-A bond, written from the
term list in README.md section 1.2 (not from ycell.py), with cross-checks.

What it does (numpy and scipy only):
  1. Builds the literal two-site bond h on C^23 x C^23 (529 x 529) and checks
     Hermiticity, exact U(1) charge, and that every move preserves the arrow
     grammar (so the penalty commutes with h).
  2. Builds the open chain [frozen M_empty] + one cell from h, on the full
     23^(n+2) space for n = 1, 2, and checks that the legal-cell subspace is
     invariant (no leakage) and equals ycell.build_cell(n).H entry by entry.
     For n = 1..4 it compares the legal-cell matrix and the port operators
     with ycell directly.
  3. Computes the second-order inter-cell hopping T in the one-particle sector
     by two new methods: (B) all orders in the conversion g, by a sparse solve
     in the two-cell charge sectors; (C) the exact g^2 coefficient
     -v_out^dag A0^-1 C A0^-1 C A0^-1 v_in, which has no scalar background.
     It compares both with ycell's one_particle_T (method A) and reports the
     eigenvalues of the Hermitian part of the phase-rotated T.
Usage: python3 bond23.py [quick]   ->  bond23_results.json
"""
import sys, os, json, itertools, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ alphabet
LET = ([('L', None, a) for a in (0, 1)] + [('R', None, a) for a in (0, 1)]
       + [('Xd', s, a) for s in 'v0+-' for a in (0, 1)]
       + [('Xx', s, a) for s in '0+-' for a in (0, 1)]
       + [('Me', None, None)] + [('M', s, None) for s in 'v0+-'])
IX = {l: i for i, l in enumerate(LET)}
D = len(LET)
assert D == 23
EXC = ("0", "+", "-")
SCH = {'v': 0, '0': 0, '+': 1, '-': -1}
ARROW = {'L': ('L', 'L'), 'R': ('R', 'R'), 'Xd': ('L', 'R'), 'Xx': ('L', 'R'),
         'Me': ('R', 'L'), 'M': ('L', 'L')}
MARKERS = [('Me', None, None)] + [('M', s, None) for s in 'v0+-']


def charge(l):
    k, s, a = l
    return (a or 0) + (SCH[s] if s is not None else 0)


def z(s, a):
    """Neutral species carries Z of the departure data bit."""
    return (-1) ** a if s == '0' else 1


WORKING = dict(t=1.0, tx=1.0, tsw=1.0, td=0.0, Delta=0.0, delta=0.0, Dx=4.0, DM=4.0,
               gp=0.05, gm=0.05, eps=0.0, musel=0.0, K=10.0, lam=1.0,
               phi=np.pi / 2,
               ph={'0': 0.0, '+': 0.0, '-': 0.0},       # stub phases
               phx={'0': 0.4, '+': 0.0, '-': 0.0},      # express phases
               phsw={'0': 0.0, '+': 0.0, '-': 0.0},     # switch phases
               phd={'0': 0.9, '+': 0.0, '-': 0.0})      # loop-link phases


def params(**kw):
    P = {k: (dict(v) if isinstance(v, dict) else v) for k, v in WORKING.items()}
    for k, v in kw.items():
        P[k] = dict(v) if isinstance(v, dict) else v
    return P


# ------------------------------------------------------------------ the bond
def two_site_kinetic(P):
    """All two-site moves (K1-K7 of README section 1.2), with h.c."""
    H = np.zeros((D * D, D * D), complex)

    def add(lo, ro, li, ri, c):          # c |lo ro><li ri| + h.c.
        o, i = IX[lo] * D + IX[ro], IX[li] * D + IX[ri]
        assert o != i
        H[o, i] += c
        H[i, o] += np.conj(c)
    for a in (0, 1):
        for b in (0, 1):
            # 1. vacuum token on the data track
            add(('L', None, a), ('Xd', 'v', b), ('Xd', 'v', a), ('R', None, b), -P['t'])
            for s in EXC:
                # 2. excited token on the data track (stub)
                add(('L', None, a), ('Xd', s, b), ('Xd', s, a), ('R', None, b),
                    -P['t'] * np.exp(1j * P['ph'][s]) * z(s, a))
                # 3. express track
                add(('L', None, a), ('Xx', s, b), ('Xx', s, a), ('R', None, b),
                    -P['tx'] * np.exp(1j * P['phx'][s]) * z(s, a))
        # 5. vacuum token into the marker
        add(('L', None, a), ('M', 'v', None), ('Xd', 'v', a), ('Me', None, None), -P['t'])
        for s in EXC:
            # 4. express track into the marker
            add(('L', None, a), ('M', s, None), ('Xx', s, a), ('Me', None, None),
                -P['tx'] * np.exp(1j * P['phx'][s]) * z(s, a))
            # 6. weak loop link, stub into the marker
            if P['td'] != 0:
                add(('L', None, a), ('M', s, None), ('Xd', s, a), ('Me', None, None),
                    -P['td'] * np.exp(1j * P['phd'][s]) * z(s, a))
            # 7. junction switch, any marker letter on the left
            for m in MARKERS:
                add(m, ('Xx', s, a), m, ('Xd', s, a),
                    -P['tsw'] * np.exp(1j * P['phsw'][s]))
    return H


def port_ops():
    """F_s on marker letters (FM) and on data-track letters (FX), 23 x 23."""
    FM, FX = {}, {}
    pairs = {'0': [('0', 'v'), ('v', '0')], '+': [('+', 'v'), ('v', '-')]}
    for f in ('0', '+'):
        A = np.zeros((D, D)); B = np.zeros((D, D))
        for (so, si) in pairs[f]:
            A[IX[('M', so, None)], IX[('M', si, None)]] = 1
            for a in (0, 1):
                B[IX[('Xd', so, a)], IX[('Xd', si, a)]] = 1
        FM[f], FX[f] = A, B
    FM['-'], FX['-'] = FM['+'].T.copy(), FX['+'].T.copy()
    return FM, FX


def port(P):
    """lambda [F0^M F0^X + e^{i phi} F+^M F-^X + e^{-i phi} F-^M F+^X]."""
    FM, FX = port_ops()
    V = (np.kron(FM['0'], FX['0']) + np.exp(1j * P['phi']) * np.kron(FM['+'], FX['-'])
         + np.exp(-1j * P['phi']) * np.kron(FM['-'], FX['+']))
    return P['lam'] * V


def onsite(P):
    """One-site operator o (README section 1.2)."""
    o = np.zeros((D, D), complex)
    dl = {'0': 0.0, '+': P['delta'], '-': -P['delta']}
    for i, (k, s, a) in enumerate(LET):
        e = 0.0
        if a is not None:
            e += P['eps'] * a
        if s in EXC:
            e += P['Delta'] + dl[s] + (P['Dx'] if k == 'Xx' else 0) + (P['DM'] if k == 'M' else 0)
        if k in ('Me', 'M'):
            e += P['musel']
        o[i, i] = e
    for (lo, li, g) in ((('Xd', '+', 0), ('Xd', '0', 1), P['gp']),
                        (('Xd', '-', 1), ('Xd', '0', 0), P['gm'])):
        o[IX[lo], IX[li]] += g
        o[IX[li], IX[lo]] += np.conj(g)
    return o


def penalty(P):
    p = np.zeros(D * D)
    for i, l in enumerate(LET):
        for j, r in enumerate(LET):
            if ARROW[l[0]][1] != ARROW[r[0]][0]:
                p[i * D + j] = P['K']
    return np.diag(p)


def bond(P):
    """The complete repeated bond h = h_kin + h_pen + lambda V + (o x 1 + 1 x o)/2."""
    o = onsite(P)
    I = np.eye(D)
    return two_site_kinetic(P) + penalty(P) + port(P) + 0.5 * (np.kron(o, I) + np.kron(I, o))


def bond_checks(P):
    h = bond(P)
    q = np.array([charge(l) for l in LET])
    Q2 = np.add.outer(q, q).reshape(-1)
    pen = np.diag(penalty(P)).real
    return dict(dim=h.shape[0],
                herm=float(np.abs(h - h.conj().T).max()),
                charge_comm=float(np.abs(h * (Q2[:, None] - Q2[None, :])).max()),
                penalty_comm=float(np.abs(h * (pen[:, None] - pen[None, :])).max()),
                nonzero=int(np.sum(np.abs(h) > 1e-15)))


# ------------------------------------------------------------------ one cell
def cell_words(n):
    """Legal cell words (sites 0..n-1 then the marker), in an independent order."""
    W = []
    for bits in itertools.product((0, 1), repeat=n):
        for j in range(n):
            for kind, sp_ in (('Xd', 'v0+-'), ('Xx', EXC)):
                for s in sp_:
                    W.append(tuple([('L', None, bits[i]) for i in range(j)] + [(kind, s, bits[j])]
                                   + [('R', None, bits[i]) for i in range(j + 1, n)]
                                   + [('Me', None, None)]))
        for s in 'v0+-':
            W.append(tuple([('L', None, b) for b in bits] + [('M', s, None)]))
    return W


def cell_from_bond(n, P):
    """H0 of one cell from the literal bond: on-site o on the n+1 cell sites, the
    kinetic two-site part on every adjacent pair including (frozen left M_empty,
    site 0).  The port is not part of H0.  Asserts closure of the legal space."""
    W = cell_words(n)
    ix = {w: i for i, w in enumerate(W)}
    hk = two_site_kinetic(P) + penalty(P)
    o = onsite(P)
    H = np.zeros((len(W), len(W)), complex)
    left = ('Me', None, None)
    for i, w in enumerate(W):
        full = (left,) + w
        for p, l in enumerate(w):
            col = o[:, IX[l]]
            for j in np.nonzero(col)[0]:
                w2 = list(w); w2[p] = LET[j]
                H[ix[tuple(w2)], i] += col[j]
        for p in range(len(full) - 1):
            col = hk[:, IX[full[p]] * D + IX[full[p + 1]]]
            for j in np.nonzero(col)[0]:
                f2 = list(full); f2[p], f2[p + 1] = LET[j // D], LET[j % D]
                assert f2[0] == left
                H[ix[tuple(f2[1:])], i] += col[j]       # KeyError = leakage
    FMc, FXc = {}, {}
    FM, FX = port_ops()
    for f in ('0', '+', '-'):
        A = np.zeros((len(W), len(W))); B = np.zeros((len(W), len(W)))
        for i, w in enumerate(W):
            for j in np.nonzero(FM[f][:, IX[w[-1]]])[0]:
                w2 = list(w); w2[-1] = LET[j]; A[ix[tuple(w2)], i] = 1
            for j in np.nonzero(FX[f][:, IX[w[0]]])[0]:
                w2 = list(w); w2[0] = LET[j]; B[ix[tuple(w2)], i] = 1
        FMc[f], FXc[f] = A, B
    Q = np.array([sum(charge(l) for l in w) for w in W])
    exc = np.array([any(l[1] in EXC for l in w) for w in W])
    data = [tuple(l[2] for l in w if l[2] is not None) for w in W]
    return dict(n=n, W=W, ix=ix, H=H, FM=FMc, FX=FXc, Q=Q, exc=exc, data=data)


def code_vectors(c, P):
    """Vacuum-walker code: E0 and G|a> from the open-path formula (eps = musel = 0)."""
    n = c['n']
    E0 = -2 * P['t'] * np.cos(np.pi / (n + 2))
    v = np.sqrt(2 / (n + 2)) * np.sin(np.pi * np.arange(1, n + 2) / (n + 2))
    G = {}
    for bits in itertools.product((0, 1), repeat=n):
        g = np.zeros(len(c['W']), complex)
        for j in range(n):
            w = tuple([('L', None, bits[i]) for i in range(j)] + [('Xd', 'v', bits[j])]
                      + [('R', None, bits[i]) for i in range(j + 1, n)] + [('Me', None, None)])
            g[c['ix'][w]] = v[j]
        g[c['ix'][tuple([('L', None, b) for b in bits] + [('M', 'v', None)])]] = v[n]
        G[bits] = g
    res = max(np.abs(c['H'] @ g - E0 * g).max() for g in G.values())
    return E0, G, res


# ------------------------------------------------------------------ cross-check with ycell
def compare_with_ycell(n, P):
    sys.path.insert(0, os.path.join(HERE, '..', 'routeA_spatial'))
    import ycell
    Py = ycell.params(t=P['t'], tx=P['tx'], tsw=P['tsw'], Delta=P['Delta'], Dx=P['Dx'], DM=P['DM'],
                      delta=P['delta'], gp=P['gp'], gm=P['gm'], phs=P['ph'], phx=P['phx'],
                      phsw=P['phsw'], bphi=P['phi'], eps=P['eps'], musel=P['musel'],
                      td=P['td'], phd=P['phd'])
    cy = ycell.build_cell(n, Py)
    c = cell_from_bond(n, P)
    perm = []
    for wy in cy['S']:
        perm.append(c['ix'][tuple(ycell.LET[l] if ycell.LET[l][0] != 'Me' else ('Me', None, None)
                                  for l in wy)])
    perm = np.array(perm)
    dH = np.abs(c['H'][np.ix_(perm, perm)] - cy['H']).max()
    dF = max(max(np.abs(c['FM'][f][np.ix_(perm, perm)] - cy['FM'][f]).max(),
                 np.abs(c['FX'][f][np.ix_(perm, perm)] - cy['FX'][f]).max()) for f in ('0', '+', '-'))
    return dict(n=n, legal_dim=len(c['W']), expected=(7 * n + 4) * 2 ** n,
                max_dH=float(dH), max_dF=float(dF))


def full_space_check(n, P):
    """Open chain [site -1, sites 0..n] with the literal bond on the full 23^(n+2)
    space (on-site o counted fully on every site).  The subspace {M_empty at
    site -1} x {legal cell} must be invariant, and the restriction must equal
    the cell Hamiltonian plus o(M_empty)."""
    L = n + 2
    hk = sp.csr_matrix(two_site_kinetic(P) + penalty(P) + port(P))
    o = sp.csr_matrix(onsite(P))
    H = None
    for p in range(L - 1):
        term = sp.kron(sp.kron(sp.identity(D ** p), hk), sp.identity(D ** (L - p - 2)))
        H = term if H is None else H + term
    for p in range(L):
        H = H + sp.kron(sp.kron(sp.identity(D ** p), o), sp.identity(D ** (L - p - 1)))
    H = H.tocsr()
    c = cell_from_bond(n, P)
    idx = []
    for w in c['W']:
        full = (('Me', None, None),) + w
        k = 0
        for l in full:
            k = k * D + IX[l]
        idx.append(k)
    idx = np.array(idx)
    mask = np.zeros(H.shape[0], bool); mask[idx] = True
    Hc = H[:, idx]
    leak = np.abs(Hc[~mask]).max() if Hc[~mask].nnz else 0.0
    Hr = Hc[idx].toarray()
    ref = c['H'] + onsite(P)[IX[('Me', None, None)], IX[('Me', None, None)]] * np.eye(len(idx))
    return dict(n=n, full_dim=D ** L, leakage=float(leak), max_diff=float(np.abs(Hr - ref).max()))


# ------------------------------------------------------------------ second order, one-particle sector
def onebody_T(n, P, order='all'):
    """T_xy = <A_x, vac_B| H2 |vac_A, B_y> with H2 = -P V Q (H0 - 2E0)^-1 Q V P,
    lambda = 1, at t_d = 0 or not.  order='all': full in g (method B).
    order='g2': exact coefficient of g_p g_m (gp = gm = g) (method C)."""
    c = cell_from_bond(n, P)
    E0, G, res = code_vectors(c, P)
    H = c['H']
    if order == 'g2':
        c0 = cell_from_bond(n, params(**dict(P, gp=0.0, gm=0.0)))
        H0c = c0['H']
        Cc = H - H0c
        Cc = Cc / P['gp']                       # gp = gm assumed
    terms = [('0', '0', 1.0), ('+', '-', np.exp(1j * P['phi'])), ('-', '+', np.exp(-1j * P['phi']))]
    zero = (0,) * n

    def e(k):
        return tuple(1 if i == k else 0 for i in range(n))

    # cell excited subspaces by charge
    sub = {q: np.where(c['exc'] & (c['Q'] == q))[0] for q in (0, 1)}
    T = np.zeros((n, n), complex)
    for (qA, qB) in ((0, 1), (1, 0)):
        iA, iB = sub[qA], sub[qB]
        dA, dB = len(iA), len(iB)
        if order == 'g2':
            HA, HB = sp.csr_matrix(H0c[np.ix_(iA, iA)]), sp.csr_matrix(H0c[np.ix_(iB, iB)])
            CA, CB = sp.csr_matrix(Cc[np.ix_(iA, iA)]), sp.csr_matrix(Cc[np.ix_(iB, iB)])
            Cpair = (sp.kron(CA, sp.identity(dB)) + sp.kron(sp.identity(dA), CB)).tocsc()
        else:
            HA, HB = sp.csr_matrix(H[np.ix_(iA, iA)]), sp.csr_matrix(H[np.ix_(iB, iB)])
        Hpair = (sp.kron(HA, sp.identity(dB)) + sp.kron(sp.identity(dA), HB)
                 - 2 * E0 * sp.identity(dA * dB)).tocsc()
        lu = spla.splu(Hpair)

        def src(a_bits, b_bits):
            v = np.zeros(dA * dB, complex)
            for (fm, fx, cf) in terms:
                va = (c['FM'][fm] @ G[a_bits])[iA]
                vb = (c['FX'][fx] @ G[b_bits])[iB]
                v += cf * np.kron(va, vb)
            return v
        vin = [src(zero, e(y)) for y in range(n)]
        vout = [src(e(x), zero) for x in range(n)]
        for y in range(n):
            if not np.any(vin[y]):
                continue
            r = lu.solve(vin[y])
            if order == 'g2':
                r = lu.solve(Cpair @ r)
                r = lu.solve(Cpair @ r)
            for x in range(n):
                T[x, y] += -np.vdot(vout[x], r)
    return T, E0, res


def herm_eigs(T):
    R = T * np.exp(-1j * np.angle(T[0, 0]))
    sc = np.abs(T).max()
    return np.linalg.eigvalsh((R + R.conj().T) / 2) / sc, float(np.abs(R.imag).max() / sc), sc


def ycell_T(n, P):
    sys.path.insert(0, os.path.join(HERE, '..', 'routeA_spatial'))
    import ycell
    Py = ycell.params(t=P['t'], tx=P['tx'], tsw=P['tsw'], Delta=P['Delta'], Dx=P['Dx'], DM=P['DM'],
                      delta=P['delta'], gp=P['gp'], gm=P['gm'], phs=P['ph'], phx=P['phx'],
                      phsw=P['phsw'], bphi=P['phi'], eps=P['eps'], musel=P['musel'],
                      td=P['td'], phd=P['phd'])
    cy = ycell.build_cell(n, Py)
    E0, G, gap = ycell.code(cy)
    H2, k = ycell.second_order(cy, G, E0)
    return ycell.one_particle_T(H2, n), gap, k


if __name__ == '__main__':
    quick = len(sys.argv) > 1 and sys.argv[1] == 'quick'
    out = {}
    P = params()
    out['bond_checks_working_point'] = bond_checks(P)
    out['bond_checks_generic'] = bond_checks(params(td=0.01, delta=0.13, eps=0.07, musel=0.3,
                                                    ph={'0': 0.2, '+': -0.3, '-': 0.5},
                                                    phsw={'0': 0.1, '+': 0.2, '-': -0.4}))
    print('bond', out['bond_checks_working_point'], out['bond_checks_generic'], flush=True)
    Pg = params(td=0.01, delta=0.13, eps=0.07, musel=0.0, gp=0.3, gm=0.2,
                ph={'0': 0.2, '+': -0.3, '-': 0.5}, phsw={'0': 0.1, '+': 0.2, '-': -0.4})
    out['cell_vs_ycell'] = [compare_with_ycell(n, Pg) for n in (1, 2, 3, 4)]
    print('cell vs ycell', out['cell_vs_ycell'], flush=True)
    out['full_space'] = [full_space_check(n, Pg) for n in ((1,) if quick else (1, 2))]
    print('full space', out['full_space'], flush=True)
    recs = []
    for (n, Dx) in ((1, 4.), (2, 4.), (3, 4.), (4, 4.), (4, 10.), (4, 20.)):
        Pn = params(Dx=Dx, DM=Dx)
        t0 = time.time()
        TA, gap, k = ycell_T(n, Pn)
        TB, E0, res = onebody_T(n, Pn, 'all')
        TC, _, _ = onebody_T(n, Pn, 'g2')
        TC = TC * Pn['gp'] ** 2
        eA, imA, sA = herm_eigs(TA)
        eB, imB, sB = herm_eigs(TB)
        eC, imC, sC = herm_eigs(TC)
        rec = dict(n=n, Dx=Dx, cell_gap=gap, poles=k, code_residual=res,
                   T_scale_A=sA, T_scale_B=sB, T_scale_g2=sC,
                   relerr_A_vs_B=float(np.abs(TA - TB).max() / sB),
                   relerr_B_vs_g2=float(np.abs(TB - TC).max() / sC),
                   herm_eigs_A=list(eA), herm_eigs_B=list(eB), herm_eigs_g2=list(eC),
                   imag_after_phase_A=imA, imag_after_phase_B=imB, imag_after_phase_g2=imC,
                   T_g2=[[str(v) for v in r] for r in TC], seconds=time.time() - t0)
        print(json.dumps({k_: rec[k_] for k_ in rec if k_ != 'T_g2'}), flush=True)
        recs.append(rec)
    out['T_working_point'] = recs
    # dependence on the express phase: T should be proportional to sin(n kappa_x)
    ph = []
    for n in (2, 3):
        for k0 in (0.4, np.pi / (2 * n)):
            Pn = params(phx={'0': k0, '+': 0.0, '-': 0.0})
            TC, _, _ = onebody_T(n, Pn, 'g2')
            ph.append(dict(n=n, phx0=k0, T_scale_g2=float(np.abs(TC).max()),
                           sin_n_kappa=float(abs(np.sin(n * k0)))))
    out['express_phase'] = ph
    print('express phase', ph, flush=True)
    with open(os.path.join(HERE, 'bond23_results.json'), 'w') as f:
        json.dump(out, f, indent=1, default=float)

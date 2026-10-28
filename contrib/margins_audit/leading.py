"""Leading (order lambda^2 g^2) logical chain of the 23-letter folded-bus cell,
in high precision, built from the cell's structure without any 2^n space.

Derivation (README section 2).  In the excited two-cell sector each cell holds
one token of species s in {0,+,-} on 2n+1 vertices d(0..n-1) (stub),
x(0..n-1) (express), M (marker).  After the diagonal +-1 gauge
|p,0,a> -> sigma(p,a)|p,0,a>, sigma = product of (-1)^{a_j} over data sites
left of p (all of them at M), the token Hamiltonian A_s no longer depends on
the data, and a conversion at stub site j becomes a Jordan-Wigner c_j or
c_j^dag of that cell times a projector on d(j).  The port creates the source
alpha^2 [P(a_L)|M,0;d0,0> + e^{i phi}|M,+;d0,-> + e^{-i phi}|M,-;d0,+>].
With H_0 - 2E_0 = A_sL x 1 + 1 x A_sR + g C, the lambda^2 g^2 coefficient is

    K = - Z^dag R0 C R0 C R0 Z ,   R0 = (A_sL x 1 + 1 x A_sR)^-1 ,

with alpha^4 and g^2 factored out.  T_xy = <A_x|K|B_y>, mu = (one-particle
diagonal blocks of both cells) - (vacuum value).  The two-carrier resolvents
are applied in the product eigenbasis of the A_s, where every restricted
source is rank one, so the cost is polynomial in n.
All literal phases of the bond are kept in A_s (no phase gauge is used).
"""
import itertools
from decimal import Decimal
import numpy as np
from hpmat import CM, ZERO, ONE, dzeros, dexpi, dcos, dsin, dpi, eigh_hp

SPEC = ['0', '+', '-']           # 0 = neutral, 1 = '+', 2 = '-'


def default_params(n):
    """C33 / C37 family: stub energy a = 4 above the vacuum (Delta = 4 - 2cos(pi/(n+2))),
    D_x = D_M = 1000, phi^x_0 = pi/(2n), port phase pi/2, C33 loop phases."""
    return dict(a=Decimal(4), Dx=Decimal(1000), DM=Decimal(1000), delta=ZERO,
                t=ONE, tx=ONE, tsw=ONE, td=ZERO,
                phs={'0': '0', '+': '0', '-': '0'},
                phx={'0': ('pi/2n', n), '+': '0', '-': '0'},
                phsw={'0': '0', '+': '0', '-': '0'},
                phd={'0': 'pi', '+': 'pi/2', '-': 'pi/2'},
                bphi='pi/2')


def phase(p):
    """phase spec -> (cos, sin) Decimals"""
    if isinstance(p, tuple) and p[0] == 'pi/2n':
        x = dpi() / (2 * p[1])
        return dcos(x), dsin(x)
    return dexpi(p)


def beta(n, a, Dx):
    """beta = 1/det(C), C the (n+1)-vertex express block, diagonal a+Dx, hopping -1."""
    d0, d1 = ONE, a + Dx
    for _ in range(n):
        d0, d1 = d1, (a + Dx) * d1 - d0
    return ONE / d1


def carrier(n, P, s):
    """Literal one-token Hamiltonian minus E_0 for species s, as a CM.
    Vertex order: d(0..n-1) = 0..n-1, x(0..n-1) = n..2n-1, M = 2n."""
    m = 2 * n + 1
    A = CM.zeros((m, m))
    ds = {'0': ZERO, '+': P['delta'], '-': -P['delta']}[s]
    for j in range(n):
        A.r[j, j] = P['a'] + ds
        A.r[n + j, n + j] = P['a'] + ds + P['Dx']
    A.r[2 * n, 2 * n] = P['a'] + ds + P['DM']

    def hop(to, frm, amp, ph):         # <to|A|frm> = -amp e^{i ph}, plus h.c.
        c, s_ = phase(ph)
        A.r[to, frm] += -amp * c
        A.i[to, frm] += -amp * s_
        A.r[frm, to] += -amp * c
        A.i[frm, to] += amp * s_
    for j in range(n - 1):
        hop(j + 1, j, P['t'], P['phs'][s])
        hop(n + j + 1, n + j, P['tx'], P['phx'][s])
    hop(2 * n, 2 * n - 1, P['tx'], P['phx'][s])        # x(n-1) -> M
    hop(n, 0, P['tsw'], P['phsw'][s])                   # d(0) -> x(0), switch
    if P['td'] != 0:
        hop(2 * n, n - 1, P['td'], P['phd'][s])         # d(n-1) -> M, loop
    return A


# conversion rules: (species, required bit) -> new species  (g = 1)
RULES = {0: [(1, 1), (0, 2)],    # 0 -> + annihilates (bit 1 -> 0); 0 -> - creates (bit 0 -> 1)
         1: [(0, 0)],            # + -> 0 creates
         2: [(1, 0)]}            # - -> 0 annihilates


def leading(n, P, verbose=False):
    """Return dict with T, mu (CM, n x n, normalized: alpha^4 g^2 lambda^2 removed),
    the vacuum scalar, and diagnostics."""
    m = 2 * n + 1
    M = 2 * n
    # one-carrier spectra, '+' and '-' share one if identical
    A = {s: carrier(n, P, SPEC[s]) for s in range(3)}
    eig = {}
    same12 = (np.all(A[1].r == A[2].r) and np.all(A[1].i == A[2].i))
    diag = {}
    for s in range(3):
        if s == 2 and same12:
            eig[2] = eig[1]
            continue
        lam, V, info = eigh_hp(A[s])
        eig[s] = (lam, V)
        diag[SPEC[s]] = dict(iterations=info['iterations'], orth=float(info['orth']),
                             offdiag=float(info['offdiag']))
        if verbose:
            print('  eig', SPEC[s], diag[SPEC[s]], flush=True)

    def cinv(s, t):
        la, lb = eig[s][0], eig[t][0]
        return np.array([[ONE / (la[k] + lb[l]) for l in range(m)] for k in range(m)], dtype=object)
    Cinv = {}
    for s in range(3):
        for t in range(3):
            Cinv[s, t] = cinv(s, t)

    cphi, sphi = dexpi(P['bphi'])
    sources = {(0, 0): None, (1, 2): (cphi, sphi), (2, 1): (cphi, -sphi)}
    # source vectors in the product eigenbasis: xhat = conj(V_sL[M,:]) conj(V_sR[d0,:])^T / C
    xhat = {}
    for src in sources:
        VL, VR = eig[src[0]][1], eig[src[1]][1]
        u = VL[M, :].conj()
        w = VR[0, :].conj()
        outer = CM(np.outer(u.r, w.r) - np.outer(u.i, w.i), np.outer(u.r, w.i) + np.outer(u.i, w.r))
        xhat[src] = outer.rmul(Cinv[src])
    # elementary restricted vectors
    elem = {}     # mid pair -> list of (src, X, j, breq, p, q)
    for src in sources:
        xs = xhat[src]
        for X in (0, 1):
            sX = src[X]
            for j in range(n):
                for breq, snew in RULES[sX]:
                    mid = (snew, src[1]) if X == 0 else (src[0], snew)
                    if X == 0:
                        p = eig[snew][1][j, :].conj()
                        row = eig[src[0]][1][j, :]
                        q = CM(np.dot(row.r, xs.r) - np.dot(row.i, xs.i), np.dot(row.r, xs.i) + np.dot(row.i, xs.r))
                    else:
                        col = eig[src[1]][1][j, :]
                        p = CM(np.dot(xs.r, col.r) - np.dot(xs.i, col.i), np.dot(xs.r, col.i) + np.dot(xs.i, col.r))
                        q = eig[snew][1][j, :].conj()
                    elem.setdefault(mid, []).append((src, X, j, breq, p, q))
    # contractions G[mid][f, i] = sum_kl conj(O_f) O_i / C_mid
    G = {}
    for mid, lst in elem.items():
        Or = np.empty((len(lst), m * m), dtype=object)
        Oi = np.empty((len(lst), m * m), dtype=object)
        for e, (_, _, _, _, p, q) in enumerate(lst):
            Or[e] = (np.outer(p.r, q.r) - np.outer(p.i, q.i)).ravel()
            Oi[e] = (np.outer(p.r, q.i) + np.outer(p.i, q.r)).ravel()
        c = Cinv[mid].ravel()
        O = CM(Or, Oi)
        Y = O.rmul(c[None, :])
        G[mid] = O.conj() @ Y.T()
    # data bookkeeping on vacuum and one-particle states
    states = [0] + [1 << x for x in range(n)] + [1 << (n + y) for y in range(n)]
    mask = (1 << n) - 1

    def popc(v):
        return bin(v).count('1')
    groups = {}
    for si, bits in enumerate(states):
        bl = bits & mask
        for mid, lst in elem.items():
            for e, (src, X, j, breq, p, q) in enumerate(lst):
                bX = (bits >> (X * n)) & mask
                if ((bX >> j) & 1) != breq:
                    continue
                sigma = -1 if popc(bX & ((1 << j) - 1)) % 2 else 1
                if src == (0, 0):
                    cr, ci = Decimal(sigma * (-1 if popc(bl) % 2 else 1)), ZERO
                else:
                    cr, ci = sources[src]
                    cr, ci = sigma * cr, sigma * ci
                nb = bits ^ (1 << (X * n + j))
                groups.setdefault((mid, nb), []).append((si, e, cr, ci))
    S = len(states)
    Kr, Ki = dzeros((S, S)), dzeros((S, S))
    for (mid, nb), lst in groups.items():
        Gm = G[mid]
        for (sf, ef, fr, fi) in lst:
            for (si, ei, ir, ii) in lst:
                # conj(cf) ci G[ef, ei]
                ar, ai = fr * ir + fi * ii, fr * ii - fi * ir
                gr, gi = Gm.r[ef, ei], Gm.i[ef, ei]
                Kr[sf, si] -= ar * gr - ai * gi
                Ki[sf, si] -= ar * gi + ai * gr
    K = CM(Kr, Ki)
    herm = (K - K.H()).absmax()
    kv = K.r[0, 0]
    T = CM(K.r[1:n + 1, n + 1:], K.i[1:n + 1, n + 1:])
    muL = CM(K.r[1:n + 1, 1:n + 1].copy(), K.i[1:n + 1, 1:n + 1].copy())
    muR = CM(K.r[n + 1:, n + 1:].copy(), K.i[n + 1:, n + 1:].copy())
    for x in range(n):
        muL.r[x, x] -= kv
        muR.r[x, x] -= kv
    mu = muL + muR
    return dict(T=T, mu=mu, muL=muL, muR=muR, Kvac=kv, herm_err=herm, eig=diag, K=K)


def alpha4(n):
    s = dsin(dpi() / (n + 2))
    return (Decimal(2) / (n + 2)) ** 2 * s ** 4

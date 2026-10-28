#!/usr/bin/env python3
"""N-flavour complex-fermion chains for the species-torus test (contrib/torus_routeB).

Flavours a = 2*alpha + beta, (alpha, beta) in {0,1}^2 (r = 2, N = r^2 = 4), as in
contrib/chat_followup/torus_encoding.md: f_a^dag carries the torus charge
w_a = e_{up,alpha} - e_{low,beta}.  Row charges R_alpha = sum_beta n_{alpha beta} and column
charges C_beta = sum_alpha n_{alpha beta} generate the torus (2r - 1 = 3 independent U(1)s).

Contents
  * local operators (one site = all N flavours, local dimension 2^N, Jordan-Wigner inside the site
    and across sites), and builders of two-site lattice terms from fermion monomials;
  * the free chain H0 = sum_a [-t_a (f^dag_{j,a} f_{j+1,a} + h.c.) - mu_a n_{j,a}], with
    t_a = 1/(2 sin k_a) and mu_a = -2 t_a cos k_a, so every flavour has Fermi points +-k_a and
    velocity exactly 1 (the one-loop flow is then the standard current-algebra one);
  * the continuum projection of quartic lattice monomials onto u(N) current-current couplings
    G_{ab,cd} psi^dag_{aR} psi_{bR} psi^dag_{cL} psi_{dL} (resonant terms only), with
    f_{x,a} -> e^{i k_a x} psi_{aR} + e^{-i k_a x} psi_{aL};
  * the one-loop flow dG_EF/dl = (1/4 pi v) f^{ACE} f^{BDF} G_AB G_CD for H = G_AB J^A_R J^B_L,
    J^A = psi^dag T^A psi, [T^A, T^B] = i f^{ABC} T^C (derived in README; SU(2) check:
    g J_R.J_L with J = psi^dag sigma/2 psi gives dg/dl = g^2/(2 pi v)), restricted to the
    resonant subspace after every step.

Written by Claude (Opus 5.5) for contrib/torus_routeB, 2026-10-27.
"""
import itertools
import numpy as np
from functools import reduce

kron = lambda *o: reduce(np.kron, o)
NF = 4


# ----------------------------------------------------------------------------- local operators
def site_ops(nf=NF):
    """Annihilators c[a] on one site (dimension 2^nf, flavour 0 most significant) and parity P."""
    a = np.array([[0, 1], [0, 0]], complex)        # basis (|0>, |1>): a|1> = |0>
    Z = np.diag([1.0, -1.0]).astype(complex)
    I = np.eye(2, dtype=complex)
    c = [kron(*([Z] * k + [a] + [I] * (nf - k - 1))) for k in range(nf)]
    P = kron(*([Z] * nf))
    return c, P


def block_ops(m, nf=NF):
    """Annihilators f[(x, a)] on a block of m sites with the global JW string."""
    c, P = site_ops(nf)
    I = np.eye(2 ** nf, dtype=complex)
    out = {}
    for x in range(m):
        for a in range(nf):
            out[(x, a)] = kron(*([P] * x + [c[a]] + [I] * (m - x - 1)))
    return out


def monomial_matrix(mono, f):
    """mono = ((dag, x, a), ...) applied left to right as written."""
    dim = f[(0, 0)].shape[0]
    M = np.eye(dim, dtype=complex)
    for (dag, x, a) in mono:
        M = M @ (f[(x, a)].conj().T if dag else f[(x, a)])
    return M


def hc(mono):
    return tuple((1 - d, x, a) for (d, x, a) in reversed(mono))


def terms_matrix(monos, m=2, nf=NF):
    f = block_ops(m, nf)
    dim = 2 ** (nf * m)
    H = np.zeros((dim, dim), complex)
    for cf, mono in monos:
        H += cf * monomial_matrix(mono, f)
    return H


# ----------------------------------------------------------------------------- torus charges
def flavour_labels(r=2):
    return [(al, be) for al in range(r) for be in range(r)]


def torus_charge_vectors(r=2):
    """w_a in Z^{2r}: components (up_0..up_{r-1}, low_0..low_{r-1})."""
    W = []
    for (al, be) in flavour_labels(r):
        w = np.zeros(2 * r, int); w[al] = 1; w[r + be] = -1
        W.append(w)
    return np.array(W)            # (N, 2r)


def torus_generators_site(r=2):
    """Species charges on one logical site: T_s = sum_a (w_a)_s n_a, s = 0..2r-1."""
    c, P = site_ops(r * r)
    n = [x.conj().T @ x for x in c]
    W = torus_charge_vectors(r)
    return [sum(W[a, s] * n[a] for a in range(r * r)) for s in range(2 * r)]


def check_torus(H, m, r=2, nf=NF):
    """max_s || [H, sum_x T_s(x)] || for an operator on m sites."""
    Ts = torus_generators_site(r)
    I = np.eye(2 ** nf)
    err = 0.0
    for T in Ts:
        Tm = sum(kron(*([I] * x + [T] + [I] * (m - x - 1))) for x in range(m))
        err = max(err, np.linalg.norm(H @ Tm - Tm @ H))
    return err


def total_number(m, nf=NF):
    c, P = site_ops(nf)
    n = sum(x.conj().T @ x for x in c)
    I = np.eye(2 ** nf)
    return sum(kron(*([I] * x + [n] + [I] * (m - x - 1))) for x in range(m))


# ----------------------------------------------------------------------------- free chain
def free_params(k):
    k = np.asarray(k, float)
    t = 1.0 / (2 * np.sin(k))
    mu = -2 * t * np.cos(k)
    return t, mu


def free_monos(k):
    t, mu = free_params(k)
    out = []
    for a in range(len(k)):
        out.append((-t[a], ((1, 0, a), (0, 1, a))))
        out.append((-t[a], ((1, 1, a), (0, 0, a))))
        out.append((-mu[a], ((1, 0, a), (0, 0, a))))
    return out


# ----------------------------------------------------------------------------- interactions
def interaction_monos(rng, g, classes, r=2):
    """Random interaction monomials (h.c. included) of the requested classes, scale g.
    Classes (all nearest-neighbour or on-site, all conserve total number and parity):
      'U'  on-site n_a n_b (a<b)                         flavour numbers conserved -> torus
      'V'  n_{0,a} n_{1,b}                               flavour numbers conserved -> torus
      'J'  exchange f+_{0a} f_{0b} f+_{1b} f_{1a} + h.c.  flavour numbers conserved -> torus
      'R'  rectangles f+_{11} f+_{22} f_{21} f_{12} on-site and across the bond      -> torus
      'P'  pair hopping f+_{0a} f+_{1a} f_{1b} f_{0b} + h.c. (a != b)               U(1) only
      'X'  other number-conserving flavour-changing quartics (correlated hopping
           f+_{0a} n_{0c} f_{1b}, on-site f+_a f+_b f_c f_d with torus charge != 0)  U(1) only
    Returns a dict class -> list of (coef, mono)."""
    N = r * r
    W = torus_charge_vectors(r)
    out = {}
    if 'U' in classes:
        out['U'] = [(g * rng.normal(), ((1, 0, a), (0, 0, a), (1, 0, b), (0, 0, b)))
                    for a in range(N) for b in range(a + 1, N)]
    if 'V' in classes:
        out['V'] = [(g * rng.normal(), ((1, 0, a), (0, 0, a), (1, 1, b), (0, 1, b)))
                    for a in range(N) for b in range(N)]
    if 'J' in classes:
        L = []
        for a in range(N):
            for b in range(N):
                if a == b:
                    continue
                cf = g * (rng.normal() + 1j * rng.normal()) / np.sqrt(2)
                m = ((1, 0, a), (0, 0, b), (1, 1, b), (0, 1, a))
                L += [(cf, m), (np.conj(cf), hc(m))]
        out['J'] = L
    if 'R' in classes:
        L = []
        # all torus-neutral, flavour-changing quartics: f+_a f+_b f_c f_d with w_a + w_b = w_c + w_d
        # and {a,b} != {c,d}; at r = 2 these are the rectangles {0,3} <-> {1,2}
        for (xa, xb, xc, xd) in [(0, 0, 0, 0), (0, 1, 1, 0), (0, 1, 0, 1), (0, 0, 1, 1), (1, 1, 0, 0)]:
            for (a, b, c_, d) in [(0, 3, 2, 1), (0, 3, 1, 2)]:
                if (xa, xb, xc, xd) == (0, 0, 0, 0) and (a, b, c_, d) == (0, 3, 1, 2):
                    continue
                assert np.all(W[a] + W[b] == W[c_] + W[d])
                cf = g * (rng.normal() + 1j * rng.normal()) / np.sqrt(2)
                m = ((1, xa, a), (1, xb, b), (0, xc, c_), (0, xd, d))
                L += [(cf, m), (np.conj(cf), hc(m))]
        out['R'] = L
    if 'P' in classes:
        L = []
        for a in range(N):
            for b in range(N):
                if a == b:
                    continue
                cf = g * (rng.normal() + 1j * rng.normal()) / np.sqrt(2)
                m = ((1, 0, a), (1, 1, a), (0, 1, b), (0, 0, b))
                L += [(cf, m), (np.conj(cf), hc(m))]
        out['P'] = L
    if 'X' in classes:
        L = []
        for a in range(N):
            for b in range(N):
                if a == b:
                    continue
                for c_ in range(N):
                    if c_ in (a, b):
                        continue
                    cf = g * (rng.normal() + 1j * rng.normal()) / np.sqrt(2) / 2
                    m = ((1, 0, a), (1, 0, c_), (0, 0, c_), (0, 1, b))
                    L += [(cf, m), (np.conj(cf), hc(m))]
        for (a, b, c_, d) in itertools.combinations(range(N), 4):
            pass
        for a, b in itertools.combinations(range(N), 2):
            for c_, d in itertools.combinations(range(N), 2):
                if (a, b) >= (c_, d):
                    continue
                if np.all(W[a] + W[b] == W[c_] + W[d]):
                    continue
                cf = g * (rng.normal() + 1j * rng.normal()) / np.sqrt(2)
                m = ((1, 0, a), (1, 0, b), (0, 0, d), (0, 0, c_))
                L += [(cf, m), (np.conj(cf), hc(m))]
        out['X'] = L
    return out


def pairing_monos(rng, g, N=NF):
    """Parity-only quadratic pairing Delta_ab f_{0a} f_{1b} + h.c. (breaks total number)."""
    L = []
    for a in range(N):
        for b in range(N):
            cf = g * (rng.normal() + 1j * rng.normal()) / np.sqrt(2)
            m = ((0, 0, a), (0, 1, b))
            L += [(cf, m), (np.conj(cf), hc(m))]
    return L


# ----------------------------------------------------------------------------- margins
def margins(k, nmax=4):
    """Smallest |sum_i s_i k_{a_i}| mod 2 pi over signed multisets of size 2..nmax that do not
    cancel identically (each flavour's net signed count nonzero for some flavour)."""
    k = np.asarray(k, float)
    N = len(k)
    best = {}
    for m in range(2, nmax + 1, 2):
        b = np.inf
        for coeffs in itertools.product(range(-m, m + 1), repeat=N):
            if sum(abs(x) for x in coeffs) > m or sum(abs(x) for x in coeffs) == 0:
                continue
            if sum(abs(x) for x in coeffs) % 2:
                continue
            p = float(np.dot(coeffs, k))
            d = abs((p + np.pi) % (2 * np.pi) - np.pi)
            b = min(b, d)
        best[m] = b
    return best


# ----------------------------------------------------------------------------- continuum projection
def perm_sign(order):
    s = 1
    o = list(order)
    for i in range(len(o)):
        for j in range(i + 1, len(o)):
            if o[i] > o[j]:
                s = -s
    return s


def project_uN(monos, k, tol=1e-8, report=False, keep_all=False):
    """Resonant part of quartic monomials as G[a,b,c,d] (coefficient of
    psi+_{aR} psi_{bR} psi+_{cL} psi_{dL}).  Also returns a list of resonant terms that are not of
    this R-L current-current form (chiral R^4/L^4 terms are dropped silently; others reported)."""
    k = np.asarray(k, float)
    N = len(k)
    G = np.zeros((N, N, N, N), complex)
    other = []
    for cf, mono in monos:
        if len(mono) != 4:
            continue
        for chir in itertools.product((1, -1), repeat=4):
            P = 0.0
            ph = 1.0 + 0j
            fields = []
            for (dag, x, a), s in zip(mono, chir):
                if dag:
                    P -= s * k[a]; ph *= np.exp(-1j * s * k[a] * x)
                else:
                    P += s * k[a]; ph *= np.exp(1j * s * k[a] * x)
                fields.append((dag, s, a))
            if not keep_all and abs((P + np.pi) % (2 * np.pi) - np.pi) > tol:
                continue
            if len(set(fields)) < 4:
                continue                       # psi psi of the same field
            nR = sum(1 for f in fields if f[1] == 1)
            if nR in (0, 4):
                continue                       # chiral: velocity renormalization only
            key = lambda i: {(1, 1): 0, (0, 1): 1, (1, -1): 2, (0, -1): 3}[(fields[i][0], fields[i][1])]
            kinds = sorted(key(i) for i in range(4))
            if kinds != [0, 1, 2, 3]:
                if abs((P + np.pi) % (2 * np.pi) - np.pi) < tol:
                    other.append((cf * ph, fields))
                continue
            order = sorted(range(4), key=key)
            sgn = perm_sign(order)
            a, b, c_, d = [fields[i][2] for i in order]
            G[a, b, c_, d] += cf * ph * sgn
    return (G, other) if report else G


def uN_basis(N):
    T = []
    for a in range(N):
        for b in range(a + 1, N):
            M = np.zeros((N, N), complex); M[a, b] = M[b, a] = 1 / np.sqrt(2); T.append(M)
            M = np.zeros((N, N), complex); M[a, b] = -1j / np.sqrt(2); M[b, a] = 1j / np.sqrt(2); T.append(M)
    for a in range(N):
        M = np.zeros((N, N), complex); M[a, a] = 1; T.append(M)
    T = np.array(T)
    nA = len(T)
    f = np.zeros((nA, nA, nA))
    for A in range(nA):
        for B in range(nA):
            Cm = T[A] @ T[B] - T[B] @ T[A]
            for C in range(nA):
                f[A, B, C] = (-1j * np.trace(Cm @ T[C])).real
    return T, f


def to_herm(Gc, T):
    """G_AB = sum G[a,b,c,d] T^A_ba T^B_dc."""
    Gh = np.einsum('abcd,Aba,Bdc->AB', Gc, T, T)
    return Gh


def to_cplx(Gh, T):
    return np.einsum('AB,Aab,Bcd->abcd', Gh, T, T)


def resonance_mask(k, tol=1e-8):
    """mask[a,b,c,d] = 1 if psi+_{aR} psi_{bR} psi+_{cL} psi_{dL} has zero momentum."""
    k = np.asarray(k, float)
    N = len(k)
    m = np.zeros((N, N, N, N))
    for a, b, c_, d in itertools.product(range(N), repeat=4):
        P = k[a] - k[b] - k[c_] + k[d]
        if abs((P + np.pi) % (2 * np.pi) - np.pi) < tol:
            m[a, b, c_, d] = 1
    return m


def torus_mask(r=2):
    W = torus_charge_vectors(r)
    N = r * r
    m = np.zeros((N, N, N, N))
    for a, b, c_, d in itertools.product(range(N), repeat=4):
        if np.all(W[a] - W[b] + W[c_] - W[d] == 0):
            m[a, b, c_, d] = 1
    return m


def classify(Gc, tol=1e-12):
    """Norms of the resonant couplings by type: Cartan (forward), pair tunnelling J^{ab}J^{ab},
    exchange J^{ab}J^{ba}, other."""
    N = Gc.shape[0]
    out = dict(cartan=0.0, pair=0.0, exchange=0.0, other=0.0)
    for a, b, c_, d in itertools.product(range(N), repeat=4):
        v = abs(Gc[a, b, c_, d]) ** 2
        if v < tol ** 2:
            continue
        if a == b and c_ == d:
            out['cartan'] += v
        elif a == c_ and b == d:
            out['pair'] += v
        elif a == d and b == c_:
            out['exchange'] += v
        else:
            out['other'] += v
    return {kk: np.sqrt(v) for kk, v in out.items()}


def momentum_matrix(k):
    """|k_a - k_b - k_c + k_d| mod 2 pi: momentum of psi+_{aR} psi_{bR} psi+_{cL} psi_{dL}."""
    k = np.asarray(k, float)
    P = k[:, None, None, None] - k[None, :, None, None] - k[None, None, :, None] + k[None, None, None, :]
    return np.abs((P + np.pi) % (2 * np.pi) - np.pi)


def one_loop(Gc0, mask, v=1.0, lmax=60.0, h=0.002, Gstop=1.0, T=None, f=None, Pmat=None, Lam0=1.0):
    """Integrate dG/dl = (1/4 pi v) f f G G in the Hermitian basis, projecting on `mask` (complex
    basis) after each step.  Stops when max|G| > max(Gstop, 5 max|G0|) ('runaway', strong
    coupling) or l > lmax ('bounded').  Returns (fate, l, G_final_cplx, history of (l, max|G|))."""
    N = Gc0.shape[0]
    if T is None:
        T, f = uN_basis(N)
    nA = len(T)
    kap = 1.0 / (4 * np.pi * v)
    fr = f.reshape(nA * nA, nA)
    # linear projector (Hermitian basis): G -> to_herm(to_cplx(G) * mask)
    Cm = np.einsum('Aab,Bcd->abcdAB', T, T).reshape(N ** 4, nA * nA)
    Hm = np.einsum('Aba,Bdc->ABabcd', T, T).reshape(nA * nA, N ** 4)
    Pm = (Hm @ (mask.reshape(-1)[:, None] * Cm)).real
    cache = {}

    def Pm_at(l):
        # scale-dependent detuning: a component of momentum p is dropped once p > Lam0 e^{-l}
        if Pmat is None:
            return Pm
        keep = ((Pmat <= Lam0 * np.exp(-l) + 1e-12) | (Pmat < 1e-9)) & (mask > 0)
        key = keep.tobytes()
        if key not in cache:
            cache[key] = (Hm @ (keep.reshape(-1)[:, None] * Cm)).real
        return cache[key]

    def B(G):
        X = np.tensordot(f, G, ([1], [1]))            # b f c
        t = np.tensordot(G, X, ([1], [0]))            # a f c
        t = t.transpose(0, 2, 1).reshape(nA * nA, nA)  # (a c) f
        return kap * (fr.T @ t)

    proj = lambda G, l=0.0: (Pm_at(l) @ G.ravel()).reshape(nA, nA)
    G = proj(to_herm(Gc0, T).real)
    l = 0.0
    hist = [(0.0, float(np.abs(G).max()))]
    gstop = max(Gstop, 5 * np.abs(G).max())
    while l < lmax:
        nG = max(np.abs(G).max(), 1e-12)
        if nG > gstop:
            return 'runaway', l, to_cplx(G, T), hist
        dl = min(h / (kap * nG), 0.05)
        k1 = B(G); k2 = B(G + 0.5 * dl * k1)
        G = proj(G + dl * k2, l + dl)
        l += dl
        if len(hist) < 2000 and (l - hist[-1][0]) > 0.25:
            hist.append((l, float(np.abs(G).max())))
    return 'bounded', l, to_cplx(G, T), hist


# ----------------------------------------------------------------------------- so(2N) projection (includes umklapp)
def field_vec(dag, a, N=NF):
    """Complex continuum field in canonical Majoranas: psi_a = (chi_2a + i chi_2a+1)/sqrt2."""
    v = np.zeros(2 * N, complex)
    v[2 * a] = 1 / np.sqrt(2)
    v[2 * a + 1] = (-1j if dag else 1j) / np.sqrt(2)
    return v


def bilinear_to_currents(f1, f2, N=NF):
    """f1 f2 (two complex fields) = sum_{A<B} r_AB J^{AB} + const, J^{AB} = i chi_A chi_B."""
    v1, v2 = field_vec(*f1, N=N), field_vec(*f2, N=N)
    pairs = list(itertools.combinations(range(2 * N), 2))
    return np.array([-1j * (v1[A] * v2[B] - v1[B] * v2[A]) for (A, B) in pairs])


def project_so(monos, k, N=NF):
    """All quartic terms with two fields at the right Fermi points (+k_a) and two at the left ones
    (-k_a), as G[AB, CD] J^{AB}_R J^{CD}_L in the so(2N) Majorana-current basis (oscillating
    terms included with their amplitude; their momentum is fixed by the charges, see so_momenta)."""
    k = np.asarray(k, float)
    nP = N * (2 * N - 1)
    G = np.zeros((nP, nP), complex)
    for cf, mono in monos:
        if len(mono) != 4:
            continue
        for chir in itertools.product((1, -1), repeat=4):
            if sum(1 for s in chir if s == 1) != 2:
                continue
            ph = 1.0 + 0j
            fields = []
            for (dag, x, a), s in zip(mono, chir):
                ph *= np.exp((-1j if dag else 1j) * s * k[a] * x)
                fields.append((dag, s, a))
            R = [i for i in range(4) if fields[i][1] == 1]
            L = [i for i in range(4) if fields[i][1] == -1]
            sgn = perm_sign(R + L)
            r = bilinear_to_currents((fields[R[0]][0], fields[R[0]][2]), (fields[R[1]][0], fields[R[1]][2]), N)
            l = bilinear_to_currents((fields[L[0]][0], fields[L[0]][2]), (fields[L[1]][0], fields[L[1]][2]), N)
            G += cf * ph * sgn * np.outer(r, l)
    return G


def so_charges(N=NF):
    """Charge eigenbasis of the so(2N) currents under the flavour numbers (see routeB_torus)."""
    c, P = site_ops(N)
    gam = []
    for x in c:
        gam += [x + x.conj().T, 1j * (x.conj().T - x)]       # gamma_2a = f + f^dag, gamma_2a+1 = -i(f - f^dag)
    pairs = list(itertools.combinations(range(2 * N), 2))
    J = [1j * gam[A] @ gam[B] for (A, B) in pairs]
    nops = [x.conj().T @ x for x in c]
    nrm = [np.trace(x.conj().T @ x).real for x in J]
    K = []
    for a in range(N):
        M = np.zeros((len(J), len(J)), complex)
        for i, JA in enumerate(J):
            comm = nops[a] @ JA - JA @ nops[a]
            for j, JB in enumerate(J):
                M[j, i] = np.trace(JB.conj().T @ comm) / nrm[j]
        K.append(M)
    coef = np.random.default_rng(7).normal(size=N) + 3.0
    w, U = np.linalg.eig(sum(cf * Km for cf, Km in zip(coef, K)))
    Uinv = np.linalg.inv(U)
    q = np.array([[(Uinv @ Km @ U)[m, m] for Km in K] for m in range(len(J))])
    assert np.abs(q.imag).max() < 1e-9
    return U, Uinv, np.round(q.real).astype(int)


def so_momenta(q, k):
    """|(q_R - q_L).k| mod 2 pi for J~_R^mu J~_L^nu (right movers at +k_a, left movers at -k_a)."""
    P = np.array([[(q[m] - q[n]) @ k for n in range(len(q))] for m in range(len(q))])
    return np.abs((P + np.pi) % (2 * np.pi) - np.pi)


def so_flow(G, f, U, Uinv, P, v=1.0, lmax=60.0, Lam0=1.0, h=0.002, Gstop=1.0):
    """One-loop flow dG/dl = (1/4 pi v) f f G G in the J basis, component (mu, nu) of Uinv G Uinv^T
    dropped once P[mu, nu] > Lam0 e^{-l}.  Returns (fate, l, G)."""
    kap = 1.0 / (4 * np.pi * v)
    nA = G.shape[0]
    fr = f.reshape(nA * nA, nA)

    def B(G):
        X = np.tensordot(f, G, ([1], [1]))
        t = np.tensordot(G, X, ([1], [0]))
        return kap * (fr.T @ t.transpose(0, 2, 1).reshape(nA * nA, nA))

    def proj(G, l):
        Gt = Uinv @ G @ Uinv.T
        Gt = Gt * (P <= Lam0 * np.exp(-l) + 1e-12)
        return np.real(U @ Gt @ U.T)
    G = proj(np.real(G), 0.0)
    gstop = max(Gstop, 5 * np.abs(G).max())
    l = 0.0
    while l < lmax:
        if np.abs(G).max() > gstop:
            return 'runaway', l, G
        dl = min(h / (kap * max(np.abs(G).max(), 1e-12)), 0.05)
        k1 = B(G); k2 = B(G + 0.5 * dl * k1)
        G = proj(G + dl * k2, l + dl)
        l += dl
    return 'bounded', l, G

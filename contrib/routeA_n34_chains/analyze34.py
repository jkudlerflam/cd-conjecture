"""Fermionic analysis of a logical bond K (4^n x 4^n, basis as in bonds34.py).

nbody(K)        exact vacuum-normal-ordered expansion K = sum_k Op_k (k-body parts)
free_chain(K)   T, M of the one-body part:  sum_j c_j^+ T c_{j+1} + h.c. + c_j^+ M c_j
fermi(T, M)     Fermi points, velocities and Bloch vectors of the free chain (complex k-roots
                of det h on |z| = 1, refined by Newton on the band that crosses zero)
hartree_fock    self-consistent Hartree-Fock one-body chain at the filling of the ground state
vertex          antisymmetrized two-body matrix element of K between plane-wave pairs
resonances      all four-fermion momentum mismatches among Fermi points, classified
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28.
"""
import itertools
import numpy as np


def jw(m):
    Z = np.diag([1.0, -1.0]); I = np.eye(2); sm = np.array([[0, 1.0], [0, 0]])
    out = []
    for j in range(m):
        M = np.ones((1, 1))
        for i in range(m):
            M = np.kron(M, Z if i < j else (sm if i == j else I))
        out.append(M.astype(complex))
    return out


_CACHE = {}
NBODY_RESIDUAL = 0.0


def ops(m):
    if m not in _CACHE:
        _CACHE[m] = jw(m)
    return _CACHE[m]


def creator(sub, m):
    """c^+_{a1} c^+_{a2} ... (increasing order) as a sparse matrix."""
    import scipy.sparse as sp
    cs = ops(m)
    M = sp.identity(2 ** m, dtype=complex, format='csr')
    for a in sub:
        M = M @ sp.csr_matrix(cs[a].conj().T)
    return M


def nbody(K):
    """Return list Op[k] (matrices) and coefficient dicts W[k][(A,B)] with
    K = sum_k sum_{|A|=|B|=k} W_AB c^+_A (c^+_B)^dag, c^+_A = prod_{a in A, increasing} c^+_a."""
    m = int(round(np.log2(len(K))))
    vac = np.zeros(2 ** m); vac[0] = 1
    res = K.copy()
    Ops, Ws = [], []
    for k in range(m + 1):
        subs = list(itertools.combinations(range(m), k))
        Cr = {A: creator(A, m) for A in subs}
        # c^+_A |0> is a single basis state (with a sign)
        kets = {A: Cr[A] @ vac for A in subs}
        idx = {A: int(np.argmax(np.abs(kets[A]))) for A in subs}
        sgn = {A: kets[A][idx[A]] for A in subs}
        W = {}
        Op = np.zeros_like(K)
        for A in subs:
            for B in subs:
                w = np.conj(sgn[A]) * res[idx[A], idx[B]] * sgn[B]
                if abs(w) > 0:
                    W[(A, B)] = w
                    Op += w * (Cr[A] @ Cr[B].conj().T).toarray()
        res = res - Op
        Ops.append(Op); Ws.append(W)
    # res is the part of K that changes the fermion number (zero for an exactly U(1) bond;
    # roundoff otherwise).  It is dropped; its relative size is kept in NBODY_RESIDUAL.
    global NBODY_RESIDUAL
    NBODY_RESIDUAL = float(np.linalg.norm(res, 2) / max(np.linalg.norm(K, 2), 1e-300))
    return Ops, Ws


def free_chain(K, n, W1=None):
    if W1 is None:
        W1 = nbody(K)[1][1]
    h = np.zeros((2 * n, 2 * n), complex)
    for ((a,), (b,)), w in W1.items():
        h[a, b] = w
    T = h[:n, n:]
    M = h[:n, :n] + h[n:, n:]
    return T, M, h


def hk(T, M, k):
    return M + T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k)


def fermi(T, M, nk=20001):
    """Fermi points of the free chain at zero energy.  Returns list of dicts
    (k, band index in sorted order, velocity dE/dk, Bloch vector)."""
    n = len(T)
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    E = np.array([np.linalg.eigvalsh(hk(T, M, k)) for k in ks])
    out = []
    for b in range(n):
        e = E[:, b]
        s = np.sign(e)
        idx = np.where(s != np.roll(s, -1))[0]
        for i in idx:
            k0, k1 = ks[i], ks[i] + 2 * np.pi / nk
            # bisection on the sorted eigenvalue b
            f0 = np.linalg.eigvalsh(hk(T, M, k0))[b]
            for _ in range(60):
                km = 0.5 * (k0 + k1)
                fm = np.linalg.eigvalsh(hk(T, M, km))[b]
                if np.sign(fm) == np.sign(f0):
                    k0, f0 = km, fm
                else:
                    k1 = km
            k = 0.5 * (k0 + k1)
            w, U = np.linalg.eigh(hk(T, M, k))
            u = U[:, b]
            dh = 1j * T * np.exp(1j * k) - 1j * T.conj().T * np.exp(-1j * k)
            v = np.vdot(u, dh @ u).real
            out.append(dict(k=float((k + np.pi) % (2 * np.pi) - np.pi), band=b, v=float(v), u=u))
    return out


def band_data(T, M, nk=4001):
    n = len(T)
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    E = np.array([np.linalg.eigvalsh(hk(T, M, k)) for k in ks])
    out = []
    for b in range(n):
        e = E[:, b]
        occ = ks[e < 0]
        centre = float(np.angle(np.exp(1j * occ).sum())) if 0 < len(occ) < nk else None
        out.append(dict(min=float(e.min()), max=float(e.max()), width=float(e.max() - e.min()),
                        filling=len(occ) / nk, centre=centre))
    return out


def ground_corr(T, M, rmax, nk=4096):
    """G_r[a,b] = <c^+_{j,a} c_{j+r,b}> of the free ground state (zero Fermi level)."""
    n = len(T)
    ks = 2 * np.pi * (np.arange(nk) + 0.5) / nk - np.pi
    G = np.zeros((2 * rmax + 1, n, n), complex)
    for k in ks:
        w, U = np.linalg.eigh(hk(T, M, k))
        P = U[:, w < 0] @ U[:, w < 0].conj().T      # <c_{k,b}^+ ... > projector, P[b,a] = <c_a^+ c_b>
        for r in range(-rmax, rmax + 1):
            G[r + rmax] += P.T * np.exp(1j * k * r)
    return G / nk


def cofactors(X):
    """Cofactor matrices of a stack X (N, k, k)."""
    N, k, _ = X.shape
    if k == 1:
        return np.ones_like(X)
    C = np.zeros_like(X)
    idx = np.arange(k)
    for p in range(k):
        for q in range(k):
            r = idx[idx != p]; c = idx[idx != q]
            C[:, p, q] = (-1) ** (p + q) * np.linalg.det(X[:, r][:, :, c])
    return C


def mean_field_onebody(Ws, G2):
    """Mean-field one-body matrix h_eff[a, b] (coefficient of c_a^+ c_b on the bond) of the operator
    sum_k sum_AB W_AB c^+_A (c^+_B)^dag in the Gaussian state with <c_a^+ c_b> = G2[a, b]
    (gradient of the Wick expectation, all bodies included)."""
    m = len(G2)
    h = np.zeros((m, m), complex)
    for k in range(1, len(Ws)):
        W = Ws[k]
        if not W:
            continue
        keys = list(W.keys())
        A = np.array([a for a, b in keys]); B = np.array([b for a, b in keys])
        w = np.array([W[key] for key in keys])
        X = G2[A[:, :, None], B[:, None, :]]
        C = cofactors(X)
        for p in range(k):
            for q in range(k):
                np.add.at(h, (A[:, p], B[:, q]), w * C[:, p, q])
    return h


def hartree_fock(K, n, iters=40, mix=0.5, tol=1e-12, W=None, verbose=False, nk=4096):
    """Self-consistent Hartree-Fock chain from all bodies of K (exact Wick gradient).
    Returns T_HF, M_HF and the history of max|change| (relative to max|T| of the one-body part)."""
    if W is None:
        W = nbody(K)[1]
    T0, M0, h0 = free_chain(K, n, W[1])
    T, M = T0.copy(), M0.copy()
    hist = []
    m = 2 * n
    for it in range(iters):
        G = ground_corr(T, M, 1, nk=nk)          # G[r+1][a,b] = <c^+_{j,a} c_{j+r,b}>
        G2 = np.zeros((m, m), complex)
        for i in range(m):
            for k in range(m):
                si, ai = divmod(i, n); sk, ak = divmod(k, n)
                G2[i, k] = G[sk - si + 1][ai, ak]
        h = mean_field_onebody(W, G2)
        Tn = 0.5 * (h[:n, n:] + h[n:, :n].conj().T)
        Mn = 0.5 * (h[:n, :n] + h[:n, :n].conj().T) + 0.5 * (h[n:, n:] + h[n:, n:].conj().T)
        ch = max(abs(Tn - T).max(), abs(Mn - M).max()) / max(abs(T0).max(), 1e-300)
        hist.append(float(ch))
        T = (1 - mix) * T + mix * Tn if it else Tn
        M = (1 - mix) * M + mix * Mn if it else Mn
        if verbose:
            print(it, ch)
        if ch < tol:
            break
    return T, M, hist


def vertex(Op2, n, ph1, ph2, ph3, ph4):
    """<1 2| Op2 |3 4> for plane waves phi_p(s, a) = u_p(a) e^{i k_p s} on the two cells,
    |3 4> = c^+(3) c^+(4)|0>.  Each ph = (k, u)."""
    m = 2 * n
    cs = ops(m)
    vac = np.zeros(2 ** m); vac[0] = 1

    def cdag(ph):
        k, u = ph
        out = np.zeros((2 ** m, 2 ** m), complex)
        for s in range(2):
            for a in range(n):
                out += u[a] * np.exp(1j * k * s) * cs[s * n + a].conj().T
        return out
    ket = cdag(ph3) @ cdag(ph4) @ vac
    bra = cdag(ph1) @ cdag(ph2) @ vac
    return np.vdot(bra, Op2 @ ket)


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def resonances(fps, top=12):
    """All p1 + p2 - p3 - p4 (mod 2 pi) over Fermi points with {1,2} != {3,4} as multisets,
    classified by the flavour (band) and chirality content.  Returns the closest ones."""
    P = len(fps)
    rows = []
    pairs = list(itertools.combinations_with_replacement(range(P), 2))
    seen = set()
    for (a, b) in pairs:
        for (c, d) in pairs:
            if (a, b) == (c, d):
                continue
            key = tuple(sorted([(a, b), (c, d)]))
            if key in seen:
                continue
            seen.add(key)
            mis = abs(wrap(fps[a]['k'] + fps[b]['k'] - fps[c]['k'] - fps[d]['k']))
            if a == b or c == d:
                if a == b and c == d:
                    # pair of identical points: Pauli forbids a local contact term; this is the
                    # 4k_F umklapp (a,a)->(c,c) with c the opposite point of the same band, or a
                    # pair transfer if the bands differ
                    kind = 'umklapp/pair (same point twice)'
                else:
                    kind = 'involves a repeated point'
            else:
                ba = {fps[x]['band'] for x in (a, b)}; bc = {fps[x]['band'] for x in (c, d)}
                if ba == bc and len(ba) == 1:
                    kind = 'intra-band (2k_F backscattering or umklapp)'
                elif len(ba) == 1 and len(bc) == 1:
                    kind = 'pair transfer'
                else:
                    kind = 'inter-band exchange / mixed'
            rows.append(dict(points=(a, b, c, d), mismatch=float(mis), kind=kind))
    rows.sort(key=lambda r: r['mismatch'])
    return rows[:top], rows

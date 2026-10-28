"""Chain specifications for the numerics.

spec string:  cell:n:param:lam:mu:level[:extra]
  cell   compile | spatial
  n      3 | 4
  param  compile: 'tau,phiT,xL' ; spatial: 'g,td,Dx'
  lam    weight of the third-order bond K3 (compile only; K = K2 + lam K3 after dividing by lambda^2)
  mu     logical field per mode, in units of max|T| ('auto' for spatial: 0.25 * 2 * min|eig R_s|)
  level  full  : the whole bond (all bodies)
         free  : the self-consistent Hartree-Fock one-body chain of the full bond
         xx    : n decoupled tight-binding chains with the Hartree-Fock Fermi points (centre and
                 width of each sorted band) and the mean |velocity| of each band
The bond is divided by max|T| of the one-body part of K2, so the widest logical band is of
order one.  Returns (K, info).
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28.
"""
import os, json, functools
import numpy as np
import analyze34 as A
from chain_report import with_field

HERE = os.path.dirname(os.path.abspath(__file__))


def bond_file(cell, n, param):
    p = param.split(',')
    if cell == 'compile':
        return os.path.join(HERE, 'data', 'bond_compile_n%d_tau%g_phi%g_xL%g.npz' % (n, float(p[0]), float(p[1]), float(p[2])))
    return os.path.join(HERE, 'data', 'bond_spatial_n%d_g%g_td%g_Dx%g.npz' % (n, float(p[0]), float(p[1]), float(p[2])))


def quad_bond(T, M, n):
    cs = A.ops(2 * n)
    H = np.zeros((4 ** n, 4 ** n), complex)
    for a in range(n):
        for b in range(n):
            hop = T[a, b] * cs[a].conj().T @ cs[n + b]
            H += hop + hop.conj().T
            H += 0.5 * M[a, b] * (cs[a].conj().T @ cs[b] + cs[n + a].conj().T @ cs[n + b])
    return 0.5 * (H + H.conj().T)


def gauge_rot(n, phi, U):
    """Many-body unitary on the bond for c_j -> e^{-i phi j} (U d)_j: returns X -> Y = Ug^dag X Ug."""
    import scipy.linalg as la
    cs = A.ops(n)
    Hg = la.logm(U)                      # U = exp(Hg), Hg anti-Hermitian
    Xop = sum(Hg[a, b] * cs[a].conj().T @ cs[b] for a in range(n) for b in range(n))
    G = la.expm(Xop)                     # c^+ -> c^+ U (one-body action on creators)
    NB = np.diag([bin(x).count('1') for x in range(2 ** n)]).astype(complex)
    D = np.kron(np.eye(2 ** n), la.expm(-1j * phi * NB))
    Ug = D @ np.kron(G, G)
    return Ug


def matched_xx(Th, Mh, n):
    fp = A.fermi(Th, Mh)
    Tx = np.zeros((n, n), complex); Mx = np.zeros((n, n), complex)
    bands = []
    for b in range(n):
        pts = [f for f in fp if f['band'] == b]
        if len(pts) != 2:
            # a band without Fermi points: a flat gapped level, filled or empty as in the chain
            e0 = np.linalg.eigvalsh(A.hk(Th, Mh, 0.0))[b]
            Mx[b, b] = -1.0 if e0 < 0 else 1.0
            bands.append(dict(band=b, gapped=True, filled=bool(e0 < 0)))
            continue
        k1, k2 = sorted(p['k'] for p in pts)
        mid = 0.5 * (k1 + k2)
        e_mid = np.linalg.eigvalsh(A.hk(Th, Mh, mid))[b]
        c = mid if e_mid < 0 else A.wrap(mid + np.pi)
        kF = 0.5 * ((k2 - k1) if e_mid < 0 else 2 * np.pi - (k2 - k1))
        v = np.mean([abs(p['v']) for p in pts])
        t = v / (2 * np.sin(kF))
        # eps(k) = M - 2t cos(k - c), Fermi points c +- kF, occupied around c
        Tx[b, b] = -t * np.exp(-1j * c)
        Mx[b, b] = 2 * t * np.cos(kF)
        bands.append(dict(band=b, centre=float(c), kF=float(kF), v=float(v)))
    Kx = quad_bond(Tx, Mx, n)
    return Kx - np.trace(Kx) / len(Kx) * np.eye(len(Kx)), bands


def bond_from_onebody_2body(h, W2, n):
    """Bond from a one-body matrix h (2n x 2n, coefficient of c_a^+ c_b) and two-body coefficients."""
    m = 2 * n
    cs = A.ops(m)
    K = np.zeros((4 ** n, 4 ** n), complex)
    for a in range(m):
        for b in range(m):
            if h[a, b] != 0:
                K += h[a, b] * cs[a].conj().T @ cs[b]
    for ((i, j), (k, l)), w in W2.items():
        K += w * cs[i].conj().T @ cs[j].conj().T @ cs[l] @ cs[k]
    return 0.5 * (K + K.conj().T)


def onebody_matrix(T, M, n):
    h = np.zeros((2 * n, 2 * n), complex)
    h[:n, n:] = T; h[n:, :n] = T.conj().T
    h[:n, :n] = 0.5 * M; h[n:, n:] = 0.5 * M
    return h


@functools.lru_cache(None)
def base(spec6):
    cell, n, param, lam, mu = spec6.split(':')[:5]
    n = int(n); lam = float(lam)
    z = np.load(bond_file(cell, n, param))
    K2 = z['K2']
    W = A.nbody(K2)[1]
    T2, M2, _ = A.free_chain(K2, n, W[1])
    s = np.abs(T2).max()
    K = K2 / s
    if cell == 'compile' and lam:
        K = K + lam * z['K3'] / s
    if mu == 'auto':
        T = T2 / s
        R = T * np.exp(-1j * np.angle(T[0, 0]))
        mu = 0.25 * 2 * np.abs(np.linalg.eigvalsh((R + R.conj().T) / 2)).min()
    mu = float(mu)
    K = with_field(K, n, mu)
    Ops, Ws = A.nbody(K)
    resid = A.NBODY_RESIDUAL
    K = sum(Ops[1:])                         # drop the scalar and the roundoff U(1)-violating part
    Th, Mh, hist = A.hartree_fock(K, n, iters=100, mix=0.5, tol=1e-12, W=Ws)
    info = dict(cell=cell, n=n, param=param, lam=lam, mu=mu, scale=float(s), u1_roundoff=resid,
                hf_change=hist[-1], body_norms=[float(np.linalg.norm(O, 2)) for O in Ops])
    return K, Ws, Th, Mh, info


@functools.lru_cache(None)
def proxy(spec6):
    """Velocity-equalized proxy.  Gauge c_j -> e^{-i Phi j} d_j (Phi = arg T_HF[0,0] + pi if the
    Hermitian part is then negative), rotate each cell to the eigenbasis of the Hermitian part R_s
    of e^{-i Phi} T_HF, and divide every normal-ordered monomial by the geometric mean of |eig R_s|
    over its fermion operators: one-body terms by sqrt(s_a s_b), two-body terms by
    (s_a s_b s_c s_d)^(1/4).  The one-body part is the Hartree-Fock chain (so the proxy has the same
    Hartree-Fock Fermi points, shifted by Phi), and the two-body part is normal ordered with respect
    to the proxy's own Hartree-Fock ground state.  Bodies >= 3 are dropped after Hartree-Fock."""
    K, Ws, Th, Mh, info = base(spec6)
    n = info['n']
    phi = float(np.angle(np.trace(Th)))
    R = np.exp(-1j * phi) * Th
    Rs = 0.5 * (R + R.conj().T)
    lam_, U = np.linalg.eigh(Rs)
    if np.all(lam_ < 0):
        phi = phi + np.pi; R = -R; Rs = -Rs; lam_, U = np.linalg.eigh(Rs)
    sc = np.abs(lam_)
    # one-body (Hartree-Fock) in the rotated, gauged basis
    Tg = U.conj().T @ R @ U
    Mg = U.conj().T @ Mh @ U
    Sm = np.diag(sc ** -0.5)
    Tp = Sm @ Tg @ Sm
    Mp = Sm @ Mg @ Sm
    # two-body coefficients in the rotated, gauged basis: transform the vacuum two-body operator
    Ug = gauge_rot(n, phi, U)
    Op2 = A.nbody(K)[0][2]
    Op2g = Ug.conj().T @ Op2 @ Ug
    W2 = A.nbody(Op2g)[1][2]
    fl = np.concatenate([np.arange(n), np.arange(n)])          # flavour of bond mode
    W2p = {}
    for ((i, j), (k, l)), w in W2.items():
        W2p[((i, j), (k, l))] = w * (sc[fl[i]] * sc[fl[j]] * sc[fl[k]] * sc[fl[l]]) ** -0.25
    # normal order the rescaled two-body part with respect to the proxy Hartree-Fock state
    G = A.ground_corr(Tp, Mp, 1)
    m = 2 * n
    G2 = np.zeros((m, m), complex)
    for i in range(m):
        for k in range(m):
            si, ai = divmod(i, n); sk, ak = divmod(k, n)
            G2[i, k] = G[sk - si + 1][ai, ak]
    hmf = A.mean_field_onebody([{}, {}, W2p], G2)
    hmf = 0.5 * (hmf + hmf.conj().T)
    hp = onebody_matrix(Tp, Mp, n) - hmf
    Kp = bond_from_onebody_2body(hp, W2p, n)
    Kp = Kp - np.trace(Kp) / len(Kp) * np.eye(len(Kp))
    pinfo = dict(info, proxy_phi=phi, proxy_scales=sc.tolist())
    return Kp, Tp, Mp, pinfo


def load(spec):
    parts = spec.split(':')
    level = parts[5]
    spec6 = ':'.join(parts[:5])
    if level.startswith('p'):
        Kp, Tp, Mp, info = proxy(spec6)
        info = dict(info, spec=spec, level=level)
        if level == 'pfull':
            return Kp, info
        if level == 'pfree':
            Kf = quad_bond(Tp, Mp, n=info['n'])
            return Kf - np.trace(Kf) / len(Kf) * np.eye(len(Kf)), info
        if level == 'pxx':
            Kx, bands = matched_xx(Tp, Mp, info['n'])
            info['xx_bands'] = bands
            return Kx, info
    K, Ws, Th, Mh, info = base(spec6)
    info = dict(info, spec=spec, level=level)
    n = info['n']
    if level == 'full':
        return K - np.trace(K) / len(K) * np.eye(len(K)), info
    if level == 'free':
        Kf = quad_bond(Th, Mh, n)
        return Kf - np.trace(Kf) / len(Kf) * np.eye(len(Kf)), info
    if level == 'xx':
        Kx, bands = matched_xx(Th, Mh, n)
        info['xx_bands'] = bands
        return Kx, info
    raise ValueError(level)


def gaussian_block_entropy(T, M, ells, nk=8192):
    """Exact block entanglement of the free ground state for blocks of ell cells."""
    n = len(T)
    rmax = max(ells)
    ks = 2 * np.pi * (np.arange(nk) + 0.5) / nk - np.pi
    G = np.zeros((rmax, n, n), complex)
    for k in ks:
        w, U = np.linalg.eigh(A.hk(T, M, k))
        P = (U[:, w < 0] @ U[:, w < 0].conj().T).T
        ph = np.exp(1j * k * np.arange(rmax))
        G += ph[:, None, None] * P[None]
    G /= nk
    out = []
    for l in ells:
        C = np.zeros((n * l, n * l), complex)
        for i in range(l):
            for j in range(l):
                r = j - i
                C[i * n:(i + 1) * n, j * n:(j + 1) * n] = G[r] if r >= 0 else G[-r].conj().T
        x = np.linalg.eigvalsh(0.5 * (C + C.conj().T))
        x = x[(x > 1e-14) & (x < 1 - 1e-14)]
        out.append(float(-(x * np.log(x) + (1 - x) * np.log(1 - x)).sum()))
    return out

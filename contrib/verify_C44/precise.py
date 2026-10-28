"""Scale-separated evaluation of the wrapped source bond C44 (11).

Plain eigendecomposition of H_(p,r)=(H_ch-e_p)xI+sW carries absolute
eigenvalue errors ~1e-16*||H|| while the slow denominators are ~s.
Here the charge operator is diagonalized first and the shift e_p is
removed exactly (gaps as eigenvalue differences). The two slow charge
states (ground and tunnelling partner) times the core are separated from
the fast charge states by an exact invariant-subspace (Riccati) iteration,
and the slow block is diagonalized at its own scale. The helix measure
is the one-flip measure with an exact zero ground energy, whose validity
is checked against the full 2^N helix in source_bond.py.
"""
import numpy as np
import scipy.linalg as sla
from core46 import Core, Params, SPECIES, port_terms
from wrap import Charge, compensator, plus
from source_bond import pair_sum, littleendian_perm


def heigh(M):
    """Hermitian eigensolver; the MRRR driver avoids a zheevd workspace
    failure seen with Accelerate on large complex matrices."""
    return sla.eigh(M, driver='evr')


def _riccati(H, nS, iters=5000, tol=1e-15):
    """Exact invariant-subspace split of Hermitian H into its first nS
    coordinates (slow) and the rest (fast). Returns (Es, Vs, Ef, Vf, res)."""
    A, B, C = H[:nS, :nS], H[:nS, nS:], H[nS:, nS:]
    Cinv = np.linalg.inv(C)
    X = -Cinv @ B.conj().T
    for it in range(iters):
        Xn = Cinv @ (X @ A + X @ B @ X - B.conj().T)
        dX = np.abs(Xn - X).max()
        X = Xn
        if dX <= tol * max(np.abs(X).max(), 1e-300):
            break
    res = float(np.abs(C @ X - X @ (A + B @ X) + B.conj().T).max())

    def sqrtm(Mh, inv=False):
        # Mh = I + G with G small and positive; use the binomial series when safe
        G = Mh - np.eye(Mh.shape[0])
        g = np.linalg.norm(G)
        if g < 1e-4:
            if inv:
                return np.eye(Mh.shape[0]) - G / 2 + 3 * (G @ G) / 8
            return np.eye(Mh.shape[0]) + G / 2 - (G @ G) / 8
        ev, ve = heigh(Mh)
        f = 1 / np.sqrt(ev) if inv else np.sqrt(ev)
        return ve @ np.diag(f) @ ve.conj().T
    nf = H.shape[0] - nS
    Ys = np.vstack([np.eye(nS), X])                    # slow invariant subspace
    Ns = Ys.conj().T @ Ys
    Ss = Ys @ sqrtm(Ns, inv=True)
    Yf = np.vstack([-X.conj().T, np.eye(nf)])          # orthogonal complement
    Nf = Yf.conj().T @ Yf
    Sf = Yf @ sqrtm(Nf, inv=True)
    Hs = Ss.conj().T @ H @ Ss
    Hf = Sf.conj().T @ H @ Sf
    Hs = (Hs + Hs.conj().T) / 2
    Hf = (Hf + Hf.conj().T) / 2
    return Hs, Ss, Hf, Sf, res, it + 1


class PreciseCell:
    """Charge levels split into three scales: ground (energy ~ s), tunnelling
    partner (~delta_1) and the rest (~U). Each block is diagonalized at its
    own norm after exact invariant-subspace separation."""

    def __init__(self, core, ch, s, bus, groups=(1, 1)):
        idx = core.sector(bus)
        self.idx = idx
        dc = len(idx)
        a, Phi = ch.eig
        a, Phi = a.copy(), Phi.copy()
        Phi[:, 0] *= np.sign(Phi[:, 0].sum())
        delta = a - a[0]
        delta[0] = 0.0
        self.delta = delta
        m = len(a)
        bonds = [b.toarray()[np.ix_(idx, idx)] for b in core.bonds]
        ends = (core.aL + core.bR).toarray()[np.ix_(idx, idx)]
        W = np.kron(np.eye(m), ends)
        for j, bj in enumerate(bonds):
            Dc = Phi.T @ np.diag(ch.Dj[:, j]) @ Phi
            W = W + np.kron(Dc, bj)
        H = np.kron(np.diag(delta), np.eye(dc)) + s * W
        n0, n1 = groups[0] * dc, (groups[0] + groups[1]) * dc
        # split {G0,G1} from the fast charge levels
        H01, S01, H2, S2, r1, i1 = _riccati(H, n1)
        # split G0 from G1 inside the effective {G0,G1} block
        H0, S0, H1, S1, r2, i2 = _riccati(H01, n0)
        E0, V0 = heigh(H0)
        E1, V1 = heigh(H1)
        E2, V2 = heigh(H2)
        self.E = np.concatenate([E0, E1, E2])
        self.U = np.hstack([S01 @ S0 @ V0, S01 @ S1 @ V1, S2 @ V2])
        self.Phi, self.dc, self.ch, self.core, self.bus = Phi, dc, ch, core, bus
        self.riccati_residual = max(r1, r2)
        self.iterations = (i1, i2)
        self.orth = float(np.abs(self.U.conj().T @ self.U - np.eye(len(self.E))).max())
        self.resid = float(np.abs(H @ self.U - self.U * self.E).max())

    def amplitudes(self, end, G):
        O = self.ch.OL if end == 'A' else self.ch.OR
        cvec_ch = self.Phi.T @ (O * self.ch.chi)
        out = {}
        for sg in SPECIES:
            F = self.core.source(end, self.bus, sg).toarray()
            cvec = (F @ G)[self.idx]
            out[sg] = self.U.conj().T @ np.kron(cvec_ch[:, None], cvec)
        return out


def helix_measure(N, M, u, Jh):
    T = np.zeros((N, N))
    for x in range(N - 1):
        T[x, x] += .5; T[x + 1, x + 1] += .5; T[x, x + 1] -= .5; T[x + 1, x] -= .5
    T[N - 1, N - 1] += 1.0
    h, v = np.linalg.eigh(T)
    uperp = np.array([-np.conj(u[1]), np.conj(u[0])])
    off = abs(uperp.conj() @ M @ u) ** 2
    w = np.concatenate([[abs(u.conj() @ M @ u) ** 2], off * np.abs(v[0]) ** 2])
    e = np.concatenate([[0.0], Jh * h])
    return np.sqrt(w), e


def precise_bonds(n, wvac, ratio=None, s=None, Jh=1.0, U=1.0, P=None, comp='small'):
    P = P or Params(n)
    core = Core(P)
    G = core.code()
    k = n + 2
    b = k + 1
    phi = 1.0 / b
    d = 2 ** n
    ch = {p: Charge(k, p, U, wvac) for p in (1, 2)}
    delta = min(ch[1].gap, ch[2].gap, Jh / ((b + 1) * (2 * b + 1)))
    if s is None:
        s = ratio * delta
    M_L, (u1, u2) = compensator(b, phi, ch[1].r, ch[2].r, comp)
    meas = {1: helix_measure(k + 1, M_L, u1, Jh), 2: helix_measure(k + 2, M_L, u2, Jh)}
    cells = {(p, r): PreciseCell(core, ch[p], s, r) for p in (1, 2) for r in (1, 2)}
    amps = {(p, r): (cells[(p, r)].amplitudes('B', G), cells[(p, r)].amplitudes('A', G))
            for p in (1, 2) for r in (1, 2)}
    perm = littleendian_perm(n)
    Pm = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            Pm[perm[i] * d + perm[j], i * d + j] = 1
    out = {}
    for p in (1, 2):
        for q in (1, 2):
            K = np.zeros((d * d, d * d), complex)
            for r in (1, 2):
                bb, terms = port_terms(P, r)
                AL = amps[(p, r)][0]
                AR = amps[(q, r)][1]
                K += pair_sum(AL, AR, cells[(p, r)].E, cells[(q, r)].E,
                              meas[p][0], meas[p][1], meas[q][0], meas[q][1], terms, bb, d)
            out['%d%d' % (p, q)] = Pm @ (s * K) @ Pm.T
    extra = dict(s=s, delta=delta, charge=ch, cells=cells, M_L=M_L,
                 riccati=max(c.riccati_residual for c in cells.values()))
    return out, extra


def kbar_bonds(n, wvac, U=1.0, P=None, dvals=None):
    """s/delta -> 0 limit: the adjacency-compressed core resolvent (C44 'Kbar')."""
    P = P or Params(n)
    core = Core(P)
    G = core.code()
    k = n + 2
    d = 2 ** n
    if dvals is None:
        ch = {p: Charge(k, p, U, wvac) for p in (1, 2)}
        dvals = {p: ch[p].adj for p in (1, 2)}
    perm = littleendian_perm(n)
    Pm = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            Pm[perm[i] * d + perm[j], i * d + j] = 1
    eig = {}
    for p in dvals:
        for r in (1, 2):
            idx = core.sector(r)
            Hc = (core.aL + core.bR).toarray()
            for j, bj in enumerate(core.bonds):
                Hc = Hc + dvals[p][j] * bj.toarray()
            Hr = Hc[np.ix_(idx, idx)]
            E, Uv = np.linalg.eigh(Hr)
            AL = {sg: Uv.conj().T @ (core.source('B', r, sg).toarray() @ G)[idx] for sg in SPECIES}
            AR = {sg: Uv.conj().T @ (core.source('A', r, sg).toarray() @ G)[idx] for sg in SPECIES}
            eig[(p, r)] = (E, AL, AR)
    out = {}
    one = (np.ones(1), np.zeros(1))
    for p in dvals:
        for q in dvals:
            K = np.zeros((d * d, d * d), complex)
            for r in (1, 2):
                bb, terms = port_terms(P, r)
                EL, AL, _ = eig[(p, r)]
                ER, _, AR = eig[(q, r)]
                K += pair_sum(AL, AR, EL, ER, one[0], one[1], one[0], one[1], terms, bb, d)
            out['%d%d' % (p, q)] = Pm @ K @ Pm.T
    return out

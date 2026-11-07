"""Single- and two-insertion crossed states for the checkerboard transfer (d = 2, XXZ family).

Physical insertion of a one-site operator O at the column between K_e and K_o at slice j.
In the U_sym = K_e^{1/2} K_o K_e^{1/2} picture the insertion state is
    |Y> = K_e^{1/2} Yhat K_o K_e^{1/2} |psi_0>   (|1> = sqrt(u_0) |psi_0>),
derived in PROOF.md section 2: the correlator of Yhat at column x with Yhat' at column x + 2r equals
<Theta Y| U_sym^{r-1} |Y'> up to the vacuum normalization, so |Y> carries the crossed spectral measure
of the insertion.  Yhat acts on the slice-j coordinate of the column (O or O^T; irrelevant for spans).
Charge: sigma^z keeps Q, sigma^+/- change Q by +-1 (sign depends on slice parity).
"""
import numpy as np
from cflib import Column, Transfer, bond_xxz_rot, momentum_basis, sector_spectrum
from scipy.linalg import eigh


def slice_op(col, vec, j, op):
    """apply single-slice operator op (2x2) on slice j to a full column vector"""
    t = vec.reshape((2,) * col.n)
    t = np.tensordot(op, t, axes=([1], [j]))
    t = np.moveaxis(t, 0, j)
    return t.reshape(-1)


class Sectors:
    """All (Q, n) sectors with eigen-decomposition of U_sym."""

    def __init__(self, Delta, beta, m, col=None, Qmax=None, h=None, mode='u1'):
        self.col = col if col is not None else Column(2, m)
        if h is None:
            h = bond_xxz_rot(Delta)
        self.T = Transfer(h, 2, beta, m, self.col)
        self.m, self.beta, self.Delta, self.mode = m, beta, Delta, mode
        Qs = self.col.charge(mode)
        self.Qs = Qs
        self.data = {}
        qlist = np.unique(Qs) if Qmax is None else [q for q in np.unique(Qs) if abs(q) <= Qmax]
        for q in qlist:
            idx = np.where(Qs == q)[0]
            for n in range(m):
                u, V, B = sector_spectrum(self.T, idx, n, vectors=True)
                if len(u) == 0:
                    continue
                self.data[(q, n)] = dict(idx=idx, B=B, u=u, V=V)
        self.u0 = max(d['u'][0] for d in self.data.values())
        d0 = self.data[(0, 0)]
        psi0 = np.zeros(self.col.N, complex)
        psi0[d0['idx']] = d0['B'] @ d0['V'][:, 0]
        self.psi0 = psi0
        assert abs(d0['u'][0] - self.u0) < 1e-12 * self.u0
        # base vector K_o K_e^{1/2} psi_0 on which insertions act
        self.base = self.T.Ko @ (self.T.Ke_half @ psi0)

    def eps(self, key):
        return -0.5 * np.log(np.clip(self.data[key]['u'], 1e-300, None) / self.u0)

    def insertion_state(self, ops):
        """ops: list of (j, 2x2 operator); returns full vector K_e^{1/2} prod Yhat K_o K_e^{1/2} psi_0"""
        v = self.base.copy()
        for j, op in ops:
            v = slice_op(self.col, v, j, op)
        return self.T.Ke_half @ v

    def project(self, vec, key):
        """coefficients of vec in the eigenbasis of sector key"""
        d = self.data[key]
        c = d['B'].conj().T @ vec[d['idx']]
        return d['V'].conj().T @ c


SZ = np.diag([1., -1.])
SP = np.array([[0., 1.], [0., 0.]])
SM = SP.T.copy()


def orth_basis(C, tol=1e-9):
    """orthonormal basis of the column span of C (SVD with relative cutoff)"""
    if C.shape[1] == 0:
        return C
    U, s, _ = np.linalg.svd(C, full_matrices=False)
    r = int((s > tol * s[0]).sum())
    return U[:, :r]


OPS1 = [('z', SZ), ('+', SP), ('-', SM)]


def single_insertion(S, slices=(0, 1)):
    """For every sector reached by one traceless insertion: orthonormal basis of the
    single-insertion span W (in the sector eigenbasis), overlaps f_j = |P_W e_j|^2,
    and the spectral weights mu_j (summed over the insertion basis, vacuum part removed)."""
    vecs = {}
    for a, op in OPS1:
        for j in slices:
            v = S.insertion_state([(j, op)])
            for key in S.data:
                c = S.project(v, key)
                if key == (0, 0):
                    c[0] = 0.0          # remove vacuum component (zero by symmetry here)
                if np.linalg.norm(c) > 1e-10 * np.linalg.norm(v):
                    vecs.setdefault(key, []).append(c)
    out = {}
    for key, cl in vecs.items():
        C = np.array(cl).T
        Qb = orth_basis(C)
        f = (abs(Qb) ** 2).sum(1)
        mu = (abs(C) ** 2).sum(1)
        out[key] = dict(f=f, mu=mu / mu.sum(), dim=Qb.shape[1], eps=S.eps(key))
    return out


def one_particle_levels(S, W):
    """one-particle levels: in each sector the dim W eigenstates with largest overlap f"""
    lev = []
    for key, w in W.items():
        top = np.argsort(-w['f'])[:w['dim']]
        for i in top:
            lev.append((w['eps'][i], key[0], key[1], w['f'][i]))
    return lev


def fock_counts(levels, m, Ecut, maxstates=3_000_000, qmod=None):
    """enumerate bosonic multisets of one-particle levels (eps, Q, n) with total eps <= Ecut;
    returns arrays of (E, Q, n mod m) for all multi-particle states including the vacuum"""
    lv = sorted(levels, key=lambda t: t[0])
    E_out, Q_out, n_out = [0.0], [0], [0]
    stack = [(0, 0.0, 0, 0)]   # (next level index, E, Q, n)
    while stack:
        i0, E, Q, n = stack.pop()
        for i in range(i0, len(lv)):
            e, q, k = lv[i][0], lv[i][1], lv[i][2]
            if E + e > Ecut:
                break
            qq = (Q + q) % qmod if qmod else Q + q
            E_out.append(E + e); Q_out.append(qq); n_out.append((n + k) % m)
            if len(E_out) > maxstates:
                raise RuntimeError('too many Fock states, lower Ecut')
            stack.append((i, E + e, qq, (n + k) % m))
    return np.array(E_out), np.array(Q_out), np.array(n_out)

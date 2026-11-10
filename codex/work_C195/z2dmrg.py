"""C190 (Claude Opus 5.5 agent, 2026-11-09): two-site DMRG for open qubit chains with exact Z2 (fermion
parity) symmetry, written for the full qubit chain at the Fibonacci cusp.

Conventions.
  Local basis s = 0 (spin up, Z = +1, empty) and s = 1 (spin down, Z = -1, one particle).
  Sp = |up><down| removes a particle, Sm = Sp^T adds one.  Total parity P = sum_j s_j mod 2, so (-1)^N = prod Z.
  Bond label q on bond b (between sites b-1 and b) = parity of the particle number on sites < b.  Labels are
  sorted on every bond (all even indices first), so each parity block is a contiguous slice.  Bond 0 has the
  single label 0 and bond L the single label P.  Every MPS tensor M[j][a, s, b] vanishes unless q_b = q_a ^ s.
  MPO channels carry the parity of the operators already placed to their left (start, Z, end even; Sp, Sm and
  the next-nearest-neighbour waiting channels odd), and W[w, v, s', s] vanishes unless m_v = m_w ^ s' ^ s.
  Environments Le[a_bra, w, a_ket] and Re[c_bra, u, c_ket] vanish unless q_bra = q_ket ^ m.

What is exploited.  The two heavy contractions of the two-site matvec (left environment with the wavefunction,
and the result with the right environment) and the environment updates are done block by block, which costs
about a quarter of the dense contraction.  The two-site SVD splits into two independent blocks.  The parity
sector is exact by construction, so no penalty term is used.

Convergence tools.  Davidson with the diagonal of the two-site problem as preconditioner (C177), with an
absolute residual target per sweep; White's density-matrix mixer (perturbation from the Sp, Sm, Z and
waiting channels, normalized, amplitude alpha per sweep) in the early sweeps, and plain SVD truncation once
alpha = 0; energy variance <H^2> - <H>^2 of the final state; checkpoints after every sweep.
"""
import numpy as np
import time
import json

I2 = np.eye(2)
ZZ_ = np.diag([1.0, -1.0])
SP = np.array([[0.0, 1.0], [0.0, 0.0]])
SM = SP.T.copy()
PX = np.array([[0.0, 1.0], [1.0, 0.0]])
PY = np.array([[0.0, -1j], [1j, 0.0]])
OPS = {'I': I2, 'Z': ZZ_, 'Sp': SP, 'Sm': SM}
OPAR = {'I': 0, 'Z': 0, 'Sp': 1, 'Sm': 1}


# ----------------------------------------------------------------------------------------------------------
# Hamiltonians
# ----------------------------------------------------------------------------------------------------------
def pauli_bond(Jx=0.0, Jy=0.0, Jz=0.0, D=0.0, Jxy_sym=0.0):
    """J_x XX + J_y YY + D (XY - YX) + J_z ZZ (+ Jxy_sym (XY + YX)), as a 4x4 matrix on |s1 s2>."""
    h = (Jx * np.kron(PX, PX) + Jy * np.kron(PY, PY) + Jz * np.kron(ZZ_, ZZ_)
         + D * (np.kron(PX, PY) - np.kron(PY, PX)) + Jxy_sym * (np.kron(PX, PY) + np.kron(PY, PX)))
    if np.max(np.abs(h.imag)) < 1e-15:
        h = h.real
    return h


def decompose(h2):
    """h2 = sum_{A,B} c[A][B] A (x) B over the operator basis I, Z, Sp, Sm (Hilbert-Schmidt projection)."""
    c = {}
    for a, A in OPS.items():
        for b, B in OPS.items():
            AB = np.kron(A, B)
            c[(a, b)] = np.trace(AB.conj().T @ h2) / np.trace(AB.conj().T @ AB)
    rec = sum(c[k] * np.kron(OPS[k[0]], OPS[k[1]]) for k in c)
    assert np.max(np.abs(rec - h2)) < 1e-12, 'operator basis decomposition failed'
    for (a, b), v in c.items():
        if abs(v) > 1e-14 and (OPAR[a] ^ OPAR[b]) == 1:
            raise ValueError(f'bond term {a}{b} = {v} breaks parity')
    return c


class MPO:
    """H = sum_j h2[j] (sites j, j+1) + sum_j h1[j] + optional range-2 terms t_pm[j] Sp_j Sm_{j+2} + t_mp[j]
    Sm_j Sp_{j+2} (no Jordan-Wigner string, used for the two-leg XX control).  A constant is kept separately."""

    def __init__(self, L, h2, h1, nnn=None):
        self.L = L
        h2 = list(h2) if isinstance(h2, (list, tuple)) else [h2] * (L - 1)
        h1 = [np.array(x, dtype=complex) for x in (h1 if isinstance(h1, (list, tuple)) else [h1] * L)]
        self.h2, self.h1_in = h2, [x.copy() for x in h1]
        coef = []
        const = 0.0
        for j, h in enumerate(h2):
            c = decompose(np.asarray(h, dtype=complex))
            for b in ('Z',):
                h1[j + 1] = h1[j + 1] + c[('I', b)] * OPS[b]
                h1[j] = h1[j] + c[(b, 'I')] * OPS[b]
            const += c[('I', 'I')].real
            coef.append(c)
        self.const = const
        has_nnn = nnn is not None
        chans = ['start', 'Z', 'end', 'Sp', 'Sm'] + (['Spw', 'Smw'] if has_nnn else [])
        self.par = np.array([0, 0, 0, 1, 1] + ([1, 1] if has_nnn else []))
        self.n0 = 3
        w = len(chans)
        self.w = w
        cx = any(np.max(np.abs(np.asarray(x).imag)) > 1e-15 for x in h1) or \
            any(abs(v.imag) > 1e-15 for c in coef for v in c.values())
        if has_nnn:
            cx = cx or any(abs(complex(t).imag) > 1e-15 for tt in nnn for t in tt)
        self.dtype = complex if cx else float
        Ws = []
        for j in range(L):
            W = np.zeros((w, w, 2, 2), dtype=complex)
            W[0, 0] = I2
            W[2, 2] = I2
            W[0, 2] = h1[j]
            if j < L - 1:
                W[0, 1] = ZZ_
                W[0, 3] = SP
                W[0, 4] = SM
            if j > 0:
                c = coef[j - 1]
                W[1, 2] = sum(c[('Z', b)] * OPS[b] for b in ('Z', 'Sp', 'Sm'))
                W[3, 2] = sum(c[('Sp', b)] * OPS[b] for b in ('Z', 'Sp', 'Sm'))
                W[4, 2] = sum(c[('Sm', b)] * OPS[b] for b in ('Z', 'Sp', 'Sm'))
            if has_nnn:
                if 0 < j < L - 1:
                    W[3, 5] = I2
                    W[4, 6] = I2
                if j > 1:
                    tpm, tmp = nnn[j - 2]
                    W[5, 2] = tpm * SM
                    W[6, 2] = tmp * SP
            if self.dtype is float:
                assert np.max(np.abs(W.imag)) < 1e-15
                W = W.real.copy()
            Ws.append(W)
        # parity selection rule of every MPO element
        for W in Ws:
            for a in range(w):
                for b in range(w):
                    for s1 in range(2):
                        for s2 in range(2):
                            if abs(W[a, b, s1, s2]) > 0:
                                assert self.par[b] == self.par[a] ^ s1 ^ s2, 'MPO parity labels violated'
        self.W = Ws

    def chs(self, m):
        return slice(0, self.n0) if m == 0 else slice(self.n0, self.w)


def gauge_phases(L, kind):
    """Site phases g_j = exp(-i theta_j) of the rotation U = prod_j exp(-i theta_j Z_j / 2).
    'stag': theta_j = (-1)^j pi/4 (C177 PROOF 1.2), which makes the symmetric-line bond real."""
    if kind in (None, 'none'):
        return np.ones(L, dtype=complex)
    if kind == 'stag':
        th = np.array([(np.pi / 4) * (1 if j % 2 == 0 else -1) for j in range(L)])
        return np.exp(-1j * th)
    raise ValueError(kind)


def rotate_bond(h2, gj, gk):
    """U h U^dag for U = u_j (x) u_k, u = diag(exp(-i theta/2), exp(+i theta/2)) with g = exp(-i theta)."""
    def u(g):
        th = -np.angle(g)
        return np.diag([np.exp(-1j * th / 2), np.exp(1j * th / 2)])
    U = np.kron(u(gj), u(gk))
    h = U @ np.asarray(h2, dtype=complex) @ U.conj().T
    if np.max(np.abs(h.imag)) < 1e-14:
        h = h.real
    return h


# ----------------------------------------------------------------------------------------------------------
# Block helpers
# ----------------------------------------------------------------------------------------------------------
def bsl(n0, D, q):
    return slice(0, n0) if q == 0 else slice(n0, D)


def bdim(n0, D, q):
    return n0 if q == 0 else D - n0


def entropy(p):
    p = np.asarray(p, dtype=float)
    p = p[p > 1e-30]
    p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


class MPS:
    def __init__(self, M, n0, P):
        self.M, self.n0, self.P = M, n0, P
        self.L = len(M)

    @property
    def dims(self):
        return [1] + [m.shape[2] for m in self.M]

    def check(self, tol=1e-13):
        """Verify the parity structure of every tensor."""
        L = self.L
        D = self.dims
        for j in range(L):
            M = self.M[j]
            for qa in (0, 1):
                for s in (0, 1):
                    qb = 1 - (qa ^ s)
                    blk = M[bsl(self.n0[j], D[j], qa), s, bsl(self.n0[j + 1], D[j + 1], qb)]
                    if blk.size and np.max(np.abs(blk)) > tol:
                        return False
        return True

    def copy(self):
        return MPS([m.copy() for m in self.M], list(self.n0), self.P)

    def save(self, fn, meta=None):
        d = {f'M{j}': m for j, m in enumerate(self.M)}
        d['n0'] = np.array(self.n0)
        d['P'] = np.array(self.P)
        if meta is not None:
            d['meta'] = np.array(json.dumps(meta))
        np.savez(fn, **d)

    @staticmethod
    def load(fn):
        z = np.load(fn, allow_pickle=False)
        L = len([k for k in z.files if k.startswith('M')])
        return MPS([z[f'M{j}'] for j in range(L)], [int(x) for x in z['n0']], int(z['P']))


def random_mps(L, P, D0=8, seed=0, dtype=float):
    rng = np.random.default_rng(seed)
    D = [1]
    for j in range(1, L):
        D.append(int(min(D0, 2 ** min(j, L - j))))
    D.append(1)
    n0 = [1] + [(d + 1) // 2 for d in D[1:L]] + [1 if P == 0 else 0]
    for j in range(1, L):
        if D[j] == 1:
            n0[j] = 1
    M = []
    for j in range(L):
        T = np.zeros((D[j], 2, D[j + 1]), dtype=dtype)
        for qa in (0, 1):
            for s in (0, 1):
                qb = qa ^ s
                a, b = bsl(n0[j], D[j], qa), bsl(n0[j + 1], D[j + 1], qb)
                shp = T[a, s, b].shape
                r = rng.standard_normal(shp)
                if dtype is complex:
                    r = r + 1j * rng.standard_normal(shp)
                T[a, s, b] = r
        M.append(T)
    psi = MPS(M, n0, P)
    right_canonicalize(psi)
    return psi


def product_mps(L, occ, dtype=float):
    occ = set(occ)
    P = len(occ) % 2
    M, n0 = [], [1]
    q = 0
    for j in range(L):
        s = 1 if j in occ else 0
        T = np.zeros((1, 2, 1), dtype=dtype)
        T[0, s, 0] = 1.0
        M.append(T)
        q ^= s
        n0.append(1 if q == 0 else 0)
    return MPS(M, n0, P)


def right_canonicalize(psi):
    """Block LQ from the right; leaves the norm in M[0] and normalizes it."""
    L = psi.L
    for j in range(L - 1, 0, -1):
        _rq_move_left(psi, j)
    psi.M[0] /= np.linalg.norm(psi.M[0])


def _rq_move_left(psi, j):
    """Make M[j] right-orthonormal (block LQ), push the R factor into M[j-1]."""
    M = psi.M[j]
    Dl, _, Dr = M.shape
    n0l, n0r = psi.n0[j], psi.n0[j + 1]
    newQ, newR, newdims = {}, {}, {}
    for ql in (0, 1):
        rows = bsl(n0l, Dl, ql)
        nr = bdim(n0l, Dl, ql)
        c0, c1 = bsl(n0r, Dr, ql), bsl(n0r, Dr, 1 - ql)   # s = 0 -> qr = ql ; s = 1 -> qr = 1 - ql
        A = np.hstack([M[rows, 0, c0], M[rows, 1, c1]])
        if nr == 0 or A.shape[1] == 0:
            newQ[ql] = np.zeros((0, A.shape[1]), dtype=M.dtype)
            newR[ql] = np.zeros((nr, 0), dtype=M.dtype)
            newdims[ql] = 0
            continue
        Q, R = np.linalg.qr(A.conj().T)          # A^dag = Q R  ->  A = R^dag Q^dag
        newQ[ql] = Q.conj().T                    # k x cols
        newR[ql] = R.conj().T                    # nr x k
        newdims[ql] = Q.shape[1]
    k0, k1 = newdims[0], newdims[1]
    Dn = k0 + k1
    Mn = np.zeros((Dn, 2, Dr), dtype=M.dtype)
    for ql in (0, 1):
        rr = slice(0, k0) if ql == 0 else slice(k0, Dn)
        c0, c1 = bsl(n0r, Dr, ql), bsl(n0r, Dr, 1 - ql)
        n_c0 = bdim(n0r, Dr, ql)
        Q = newQ[ql]
        Mn[rr, 0, c0] = Q[:, :n_c0]
        Mn[rr, 1, c1] = Q[:, n_c0:]
    Rfull = np.zeros((Dl, Dn), dtype=M.dtype)
    Rfull[bsl(n0l, Dl, 0), 0:k0] = newR[0]
    Rfull[bsl(n0l, Dl, 1), k0:Dn] = newR[1]
    psi.M[j] = Mn
    psi.M[j - 1] = np.tensordot(psi.M[j - 1], Rfull, axes=([2], [0]))
    psi.n0[j] = k0


def _qr_move_right(psi, j):
    """Make M[j] left-orthonormal (block QR), push R into M[j+1]."""
    M = psi.M[j]
    Dl, _, Dr = M.shape
    n0l, n0r = psi.n0[j], psi.n0[j + 1]
    newQ, newR, newdims = {}, {}, {}
    for qr in (0, 1):
        r0, r1 = bsl(n0l, Dl, qr), bsl(n0l, Dl, 1 - qr)   # s = 0: ql = qr ; s = 1: ql = 1 - qr
        cols = bsl(n0r, Dr, qr)
        A = np.vstack([M[r0, 0, cols], M[r1, 1, cols]])
        nc = bdim(n0r, Dr, qr)
        if nc == 0 or A.shape[0] == 0:
            newQ[qr] = np.zeros((A.shape[0], 0), dtype=M.dtype)
            newR[qr] = np.zeros((0, nc), dtype=M.dtype)
            newdims[qr] = 0
            continue
        Q, R = np.linalg.qr(A)
        newQ[qr], newR[qr], newdims[qr] = Q, R, Q.shape[1]
    k0, k1 = newdims[0], newdims[1]
    Dn = k0 + k1
    Mn = np.zeros((Dl, 2, Dn), dtype=M.dtype)
    for qr in (0, 1):
        cc = slice(0, k0) if qr == 0 else slice(k0, Dn)
        r0, r1 = bsl(n0l, Dl, qr), bsl(n0l, Dl, 1 - qr)
        n_r0 = bdim(n0l, Dl, qr)
        Mn[r0, 0, cc] = newQ[qr][:n_r0]
        Mn[r1, 1, cc] = newQ[qr][n_r0:]
    Rfull = np.zeros((Dn, Dr), dtype=M.dtype)
    Rfull[0:k0, bsl(n0r, Dr, 0)] = newR[0]
    Rfull[k0:Dn, bsl(n0r, Dr, 1)] = newR[1]
    psi.M[j] = Mn
    psi.M[j + 1] = np.tensordot(Rfull, psi.M[j + 1], axes=([1], [0]))
    psi.n0[j + 1] = k0


def flip_site(psi, j):
    """Apply X_j: changes the parity sector; labels of bonds > j are swapped (blocks reordered)."""
    L = psi.L
    D = psi.dims
    M = [m.copy() for m in psi.M]
    M[j] = M[j][:, ::-1, :].copy()
    n0 = list(psi.n0)
    for b in range(j + 1, L + 1):
        n0b, Db = n0[b], D[b]
        perm = np.r_[np.arange(n0b, Db), np.arange(0, n0b)]
        M[b - 1] = M[b - 1][:, :, perm]
        if b < L:
            M[b] = M[b][perm, :, :]
        n0[b] = Db - n0b
    return MPS(M, n0, psi.P ^ 1)


# ----------------------------------------------------------------------------------------------------------
# Environments (block-sparse contractions on dense storage)
# ----------------------------------------------------------------------------------------------------------
def env_left(Le, A, W, n0a, n0b, mpo):
    """Lnew[b', v, b] = sum conj(A[a', s', b']) Le[a', w, a] W[w, v, s', s] A[a, s, b]."""
    Da, _, Db = A.shape
    w = mpo.w
    dt = np.result_type(Le.dtype, A.dtype, W.dtype)
    # T1[a', w, s, b] = sum_a Le[a', w, a] A[a, s, b]   (block over a and m_w)
    T1 = np.zeros((Da, w, 2, Db), dtype=dt)
    for qa in (0, 1):
        a = bsl(n0a, Da, qa)
        na = bdim(n0a, Da, qa)
        if na == 0:
            continue
        cols = []
        for s in (0, 1):
            cols.append(A[a, s, bsl(n0b, Db, qa ^ s)])
        C = np.hstack(cols)
        for m in (0, 1):
            ap = bsl(n0a, Da, qa ^ m)
            nap = bdim(n0a, Da, qa ^ m)
            ch = mpo.chs(m)
            nm = ch.stop - ch.start
            if nap == 0 or nm == 0 or C.shape[1] == 0:
                continue
            Lb = Le[ap, ch, a].reshape(nap * nm, na)
            R = (Lb @ C).reshape(nap, nm, -1)
            off = 0
            for s in (0, 1):
                nb = bdim(n0b, Db, qa ^ s)
                T1[ap, ch, s, bsl(n0b, Db, qa ^ s)] = R[:, :, off:off + nb]
                off += nb
    # T2[a', b, v, s'] = sum_{w, s} T1[a', w, s, b] W[w, v, s', s]
    T2 = np.tensordot(T1, W, axes=([1, 2], [0, 3]))
    # Lnew[b', v, b] = sum_{a', s'} conj(A[a', s', b']) T2[a', b, v, s']
    Ln = np.zeros((Db, w, Db), dtype=dt)
    for qbp in (0, 1):
        bp = bsl(n0b, Db, qbp)
        nbp = bdim(n0b, Db, qbp)
        if nbp == 0:
            continue
        for m in (0, 1):
            qb = qbp ^ m
            bk = bsl(n0b, Db, qb)
            nb = bdim(n0b, Db, qb)
            ch = mpo.chs(m)
            nm = ch.stop - ch.start
            if nb == 0 or nm == 0:
                continue
            acc = np.zeros((nbp, nb * nm), dtype=dt)
            for sp in (0, 1):
                qap = qbp ^ sp
                ap = bsl(n0a, Da, qap)
                if bdim(n0a, Da, qap) == 0:
                    continue
                Ab = A[ap, sp, bp]                         # (na', nb')
                Tb = T2[ap, bk, ch, sp]                    # (na', nb, nm)
                acc += Ab.conj().T @ Tb.reshape(Tb.shape[0], -1)
            Ln[bp, ch, bk] = acc.reshape(nbp, nb, nm).transpose(0, 2, 1)
    return Ln


def env_right(Re, B, W, n0a, n0b, mpo):
    """Rnew[a', w, a] = sum conj(B[a', s', b']) W[w, v, s', s] B[a, s, b] Re[b', v, b]."""
    Da, _, Db = B.shape
    w = mpo.w
    dt = np.result_type(Re.dtype, B.dtype, W.dtype)
    # T1[a, s, b', v] = sum_b B[a, s, b] Re[b', v, b]
    T1 = np.zeros((Da, 2, Db, w), dtype=dt)
    for qb in (0, 1):
        b = bsl(n0b, Db, qb)
        nb = bdim(n0b, Db, qb)
        if nb == 0:
            continue
        rows = []
        for s in (0, 1):
            rows.append(B[bsl(n0a, Da, qb ^ s), s, b])
        Rw = np.vstack(rows)                               # (sum na, nb)
        for m in (0, 1):
            qbp = qb ^ m
            bp = bsl(n0b, Db, qbp)
            nbp = bdim(n0b, Db, qbp)
            ch = mpo.chs(m)
            nm = ch.stop - ch.start
            if nbp == 0 or nm == 0 or Rw.shape[0] == 0:
                continue
            Rb = Re[bp, ch, b].reshape(nbp * nm, nb)       # (b' v, b)
            X = (Rw @ Rb.T).reshape(-1, nbp, nm)            # (a s, b', v)
            off = 0
            for s in (0, 1):
                na = bdim(n0a, Da, qb ^ s)
                T1[bsl(n0a, Da, qb ^ s), s, bp, ch] = X[off:off + na]
                off += na
    # T2[a, b', w, s'] = sum_{s, v} T1[a, s, b', v] W[w, v, s', s]
    T2 = np.tensordot(T1, W, axes=([1, 3], [3, 1]))       # a b' w s'
    Rn = np.zeros((Da, w, Da), dtype=dt)
    for qap in (0, 1):
        ap = bsl(n0a, Da, qap)
        nap = bdim(n0a, Da, qap)
        if nap == 0:
            continue
        for m in (0, 1):
            qa = qap ^ m
            ak = bsl(n0a, Da, qa)
            na = bdim(n0a, Da, qa)
            ch = mpo.chs(m)
            nm = ch.stop - ch.start
            if na == 0 or nm == 0:
                continue
            acc = np.zeros((nap, na * nm), dtype=dt)
            for sp in (0, 1):
                qbp = qap ^ sp
                bp = bsl(n0b, Db, qbp)
                if bdim(n0b, Db, qbp) == 0:
                    continue
                Bb = B[ap, sp, bp]                         # (na', nb')
                Tb = T2[ak, bp, ch, sp]                    # (na, nb', nm)
                acc += Bb.conj() @ Tb.transpose(1, 2, 0).reshape(Tb.shape[1], -1)
            Rn[ap, ch, ak] = acc.reshape(nap, nm, na)
    return Rn


# ----------------------------------------------------------------------------------------------------------
# Two-site effective Hamiltonian
# ----------------------------------------------------------------------------------------------------------
class TwoSite:
    def __init__(self, Le, W1, W2, Re, n0a, Da, n0c, Dc, mpo, P_left_check=None):
        self.n0a, self.Da, self.n0c, self.Dc, self.mpo = n0a, Da, n0c, Dc, mpo
        self.dt = np.result_type(Le.dtype, W1.dtype, Re.dtype)
        self.blocks = []
        off = 0
        for qa in (0, 1):
            for s1 in (0, 1):
                for s2 in (0, 1):
                    qc = qa ^ s1 ^ s2
                    shp = (bdim(n0a, Da, qa), bdim(n0c, Dc, qc))
                    n = shp[0] * shp[1]
                    self.blocks.append((qa, s1, s2, qc, shp, off, n))
                    off += n
        self.size = off
        self.Lb = {}
        for qa in (0, 1):
            na = bdim(n0a, Da, qa)
            for m in (0, 1):
                ch = mpo.chs(m)
                nm = ch.stop - ch.start
                nap = bdim(n0a, Da, qa ^ m)
                if na and nap and nm:
                    self.Lb[(qa, m)] = np.ascontiguousarray(
                        Le[bsl(n0a, Da, qa ^ m), ch, bsl(n0a, Da, qa)].reshape(nap * nm, na))
        self.RbT = {}
        for qcp in (0, 1):
            ncp = bdim(n0c, Dc, qcp)
            for m in (0, 1):
                ch = mpo.chs(m)
                nm = ch.stop - ch.start
                nc = bdim(n0c, Dc, qcp ^ m)
                if ncp and nc and nm:
                    Rb = Re[bsl(n0c, Dc, qcp), ch, bsl(n0c, Dc, qcp ^ m)]   # (c', u, c)
                    self.RbT[(qcp, m)] = np.ascontiguousarray(Rb.reshape(ncp, nm * nc).T)
        self.W12 = np.einsum('wvas,vubt->wuabst', W1, W2)               # (w, u, s1', s2', s1, s2)
        wd = mpo.w
        # (u s1' s2') x (w s1 s2) matrix for the batched step 2
        self.W12m = np.ascontiguousarray(self.W12.transpose(1, 2, 3, 0, 4, 5).reshape(wd * 4, wd * 4))
        self.W1 = W1
        self.W2 = W2
        self.Le, self.Re = Le, Re
        # diagonal
        l = np.einsum('awa->aw', Le)
        r = np.einsum('cuc->cu', Re)
        w1 = np.einsum('wvss->wvs', W1)
        w2 = np.einsum('vuss->vus', W2)
        t = np.einsum('aw,wvs->avs', l, w1)
        t = np.einsum('avs,vut->astu', t, w2)
        dg = np.einsum('astu,cu->astc', t, r)
        self.diag = self.pack_dense(dg).real.copy()

    def build_kron(self):
        """Preconditioner H0 = H_{L,j} (x) 1 + 1 (x) H_{j+1,R} restricted to each parity block r, where H_{L,j} is
        the exact Hamiltonian of the left block plus site j, H_{j+1,R} that of site j+1 plus the right block,
        and the central bond contributes only its separable part c_zz (Z_j + Z_{j+1} - 1) (the remainder of the
        central ZZ is 4 c_zz n_j n_{j+1}, nonzero only on states that cost V).  The pairing and hopping across
        the central bond are left out."""
        n0a, Da, n0c, Dc = self.n0a, self.Da, self.n0c, self.Dc
        W1, W2, Le, Re = self.W1, self.W2, self.Le, self.Re
        czz = W2[1, 2, 0, 0]          # W1[0,1] = Z and W2[1,2] = c_zz Z
        HL = np.einsum('awb,wst->asbt', Le, W1[:, 2]) + czz * np.einsum('ab,st->asbt', np.eye(Da), ZZ_)
        HR = np.einsum('wst,cwd->sctd', W2[0, :], Re) + czz * np.einsum('st,cd->sctd', ZZ_, np.eye(Dc))
        HR = HR - czz * np.einsum('st,cd->sctd', np.eye(2), np.eye(Dc))
        self.kron = {}
        for r in (0, 1):
            a0, a1 = bsl(n0a, Da, r), bsl(n0a, Da, 1 - r)
            c0, c1 = bsl(n0c, Dc, r), bsl(n0c, Dc, 1 - r)
            rows = [(a0, 0), (a1, 1)]
            cols = [(0, c0), (1, c1)]
            HLr = np.block([[HL[ra, rs, ca, cs] for (ca, cs) in rows] for (ra, rs) in rows])
            HRr = np.block([[HR[rs, rc, cs, cc] for (cs, cc) in cols] for (rs, rc) in cols])
            if HLr.size == 0 or HRr.size == 0:
                self.kron[r] = None
                continue
            lam, U = np.linalg.eigh((HLr + HLr.conj().T) / 2)
            mu, V = np.linalg.eigh((HRr + HRr.conj().T) / 2)
            self.kron[r] = (lam, U, mu, V)
        self.kron_min = min(k[0][0] + k[2][0] for k in self.kron.values() if k is not None)

    def precond(self, r_vec, theta, sigma):
        th = self.unpack_dense(r_vec)
        out = np.zeros_like(th)
        n0a, Da, n0c, Dc = self.n0a, self.Da, self.n0c, self.Dc
        for r in (0, 1):
            if self.kron.get(r) is None:
                continue
            lam, U, mu, V = self.kron[r]
            B = _theta_block(th, n0a, Da, n0c, Dc, r)
            X = U.conj().T @ B @ V.conj()
            den = lam[:, None] + mu[None, :] - min(theta, self.kron_min) + sigma
            X = U @ (X / den) @ V.T
            na0, nc0 = bdim(n0a, Da, r), bdim(n0c, Dc, r)
            a0, a1 = bsl(n0a, Da, r), bsl(n0a, Da, 1 - r)
            c0, c1 = bsl(n0c, Dc, r), bsl(n0c, Dc, 1 - r)
            out[a0, 0, 0, c0] = X[:na0, :nc0]
            out[a0, 0, 1, c1] = X[:na0, nc0:]
            out[a1, 1, 0, c0] = X[na0:, :nc0]
            out[a1, 1, 1, c1] = X[na0:, nc0:]
        return self.pack_dense(out)

    def build_kron4(self):
        """Preconditioner block-diagonal in (s1, s2): on the block (qa, s1, s2, qc),
        H0 = H_L^(s1)[qa] (x) 1 + 1 (x) H_R^(s2)[qc] + c(s1, s2), with H_L^(s1) = <s1| H_{left block + site j} |s1>
        (exact block Hamiltonian, the diagonal part in s1 of every term touching site j), H_R^(s2) likewise for
        site j+1 plus the right block, and c the diagonal of the central bond (so a doubly occupied bond costs V).
        Only the terms that flip s1 or s2 (pairing and hopping on the three bonds touching the two sites) are left
        out.  Applied as (H0 - E0 + sigma)^-1 with E0 = min(theta, lowest eigenvalue of H0)."""
        n0a, Da, n0c, Dc, mpo = self.n0a, self.Da, self.n0c, self.Dc, self.mpo
        W1, W2, Le, Re = self.W1, self.W2, self.Le, self.Re
        ev = mpo.chs(0)
        self.k4L, self.k4R = {}, {}
        for s1 in (0, 1):
            for qa in (0, 1):
                a = bsl(n0a, Da, qa)
                if bdim(n0a, Da, qa) == 0:
                    continue
                H = np.tensordot(Le[a, ev, a], W1[ev, 2, s1, s1], axes=([1], [0]))
                lam, U = np.linalg.eigh((H + H.conj().T) / 2)
                self.k4L[(s1, qa)] = (lam, U)
        for s2 in (0, 1):
            for qc in (0, 1):
                c = bsl(n0c, Dc, qc)
                if bdim(n0c, Dc, qc) == 0:
                    continue
                H = np.tensordot(Re[c, ev, c], W2[0, ev, s2, s2], axes=([1], [0]))
                mu, V = np.linalg.eigh((H + H.conj().T) / 2)
                self.k4R[(s2, qc)] = (mu, V)
        self.k4c = {(s1, s2): float(np.real(W1[0, 1, s1, s1] * W2[1, 2, s2, s2]))       # central ZZ only
                    for s1 in (0, 1) for s2 in (0, 1)}
        e0 = np.inf
        for qa, s1, s2, qc, shp, off, n in self.blocks:
            if n:
                e0 = min(e0, self.k4L[(s1, qa)][0][0] + self.k4R[(s2, qc)][0][0] + self.k4c[(s1, s2)])
        self.k4e0 = e0

    def precond4(self, r_vec, theta, sigma):
        out = np.empty_like(r_vec)
        sh = min(theta, self.k4e0) - sigma
        for qa, s1, s2, qc, shp, off, n in self.blocks:
            if n == 0:
                continue
            lam, U = self.k4L[(s1, qa)]
            mu, V = self.k4R[(s2, qc)]
            X = U.conj().T @ r_vec[off:off + n].reshape(shp) @ V.conj()
            X = X / (lam[:, None] + mu[None, :] + (self.k4c[(s1, s2)] - sh))
            out[off:off + n] = (U @ X @ V.T).reshape(-1)
        return out

    def pack_dense(self, th):
        out = np.empty(self.size, dtype=th.dtype)
        for qa, s1, s2, qc, shp, off, n in self.blocks:
            if n:
                out[off:off + n] = th[bsl(self.n0a, self.Da, qa), s1, s2, bsl(self.n0c, self.Dc, qc)].reshape(-1)
        return out

    def unpack_dense(self, x):
        th = np.zeros((self.Da, 2, 2, self.Dc), dtype=x.dtype)
        for qa, s1, s2, qc, shp, off, n in self.blocks:
            if n:
                th[bsl(self.n0a, self.Da, qa), s1, s2, bsl(self.n0c, self.Dc, qc)] = x[off:off + n].reshape(shp)
        return th

    def step1(self, x):
        """Y[a', w, s1, s2, c] = sum_a Le[a', w, a] theta[a, s1, s2, c]."""
        n0a, Da, n0c, Dc, mpo = self.n0a, self.Da, self.n0c, self.Dc, self.mpo
        dt = np.result_type(self.dt, x.dtype)
        Y = np.zeros((Da, mpo.w, 2, 2, Dc), dtype=dt)
        for qa in (0, 1):
            bl = [b for b in self.blocks if b[0] == qa and b[6] > 0]
            if not bl:
                continue
            na = bl[0][4][0]
            C = np.hstack([x[b[5]:b[5] + b[6]].reshape(b[4]) for b in bl])
            for m in (0, 1):
                if (qa, m) not in self.Lb:
                    continue
                ch = mpo.chs(m)
                nm = ch.stop - ch.start
                nap = bdim(n0a, Da, qa ^ m)
                R = (self.Lb[(qa, m)] @ C).reshape(nap, nm, -1)
                off = 0
                for b in bl:
                    nc = b[4][1]
                    Y[bsl(n0a, Da, qa ^ m), ch, b[1], b[2], bsl(n0c, Dc, b[3])] = R[:, :, off:off + nc]
                    off += nc
        return Y

    def matvec(self, x):
        n0a, Da, n0c, Dc, mpo = self.n0a, self.Da, self.n0c, self.Dc, self.mpo
        Y = self.step1(x)
        wd = mpo.w
        Zt = np.matmul(self.W12m.astype(np.result_type(self.W12m, Y), copy=False), Y.reshape(Da, wd * 4, Dc))
        Zt = Zt.reshape(Da, wd, 2, 2, Dc)                                 # (a', u, s1', s2', c)
        out = np.empty(self.size, dtype=Zt.dtype)
        for qap, s1, s2, qcp, shp, off, n in self.blocks:
            if n == 0:
                continue
            acc = np.zeros(shp, dtype=Zt.dtype)
            for m in (0, 1):
                if (qcp, m) not in self.RbT:
                    continue
                ch = mpo.chs(m)
                Zb = Zt[bsl(n0a, Da, qap), ch, s1, s2, bsl(n0c, Dc, qcp ^ m)]   # (a', u, c)
                Zb = Zb.reshape(shp[0], Zb.shape[1] * Zb.shape[2])
                acc += Zb @ self.RbT[(qcp, m)]
            out[off:off + n] = acc.reshape(-1)
        return out

    def mixer_rho_left(self, x, r):
        """Perturbation for the left density matrix of block r: X X^dag with X[(a', s1'), (s2, c, v)] =
        (Le W1 theta) over the non-identity channels v (Z, Sp, Sm, waiting)."""
        n0a, Da = self.n0a, self.Da
        Y = self.step1(x)
        Xt = np.tensordot(Y, self.W1, axes=([1, 2], [0, 3]))           # (a', s2, c, v, s1')
        keep = [v for v in range(self.mpo.w) if v not in (0, 2)]
        Xt = Xt[:, :, :, keep, :]
        ncol = int(np.prod(Xt.shape[1:4]))
        rows = [Xt[bsl(n0a, Da, r), :, :, :, 0].reshape(bdim(n0a, Da, r), ncol),
                Xt[bsl(n0a, Da, 1 - r), :, :, :, 1].reshape(bdim(n0a, Da, 1 - r), ncol)]
        X = np.vstack(rows)
        return X @ X.conj().T

    def mixer_rho_right(self, x, r):
        """Same for the right density matrix of block r: columns (s2', c'), X = theta W2 Re over channels v."""
        n0c, Dc = self.n0c, self.Dc
        th = self.unpack_dense(x)
        T = np.tensordot(th, self.Re, axes=([3], [2]))                  # (a, s1, s2, c', u)
        T = np.tensordot(T, self.W2, axes=([2, 4], [3, 1]))            # (a, s1, c', v, s2')
        keep = [v for v in range(self.mpo.w) if v not in (0, 2)]
        T = T[:, :, :, keep, :]
        cols = [T[:, :, bsl(n0c, Dc, r), :, 0], T[:, :, bsl(n0c, Dc, 1 - r), :, 1]]   # s2' = 0: qc' = r
        cols = [c.transpose(2, 0, 1, 3).reshape(c.shape[2], c.shape[0] * c.shape[1] * c.shape[3]) for c in cols]
        X = np.vstack(cols)                                              # (cols of block r, rest)
        # same convention as Bm^dag Bm in split(): conj of the right reduced density matrix
        return X.conj() @ X.T


def davidson(mv, x0, diag, res_tol=1e-8, maxit=400, mmax=32, keep=4, sigma=1e-4, prec=None):
    n = x0.size
    dt = np.result_type(x0.dtype, float)
    x = x0 / np.linalg.norm(x0)
    mmax = min(mmax, n)
    keep = min(keep, mmax - 1)
    V = np.zeros((mmax, n), dtype=dt)
    AV = np.zeros((mmax, n), dtype=dt)
    V[0] = x
    AV[0] = mv(x)
    m = 1
    nmv = 1
    theta, u, rn = None, x, np.inf
    for it in range(maxit):
        Hs = V[:m].conj() @ AV[:m].T
        Hs = (Hs + Hs.conj().T) / 2
        th, Y = np.linalg.eigh(Hs)
        y = Y[:, 0]
        theta = th[0]
        u = y @ V[:m]
        Au = y @ AV[:m]
        r = Au - theta * u
        rn = np.linalg.norm(r)
        if rn < res_tol or m >= n:
            break
        if prec is not None:
            t = prec(r, theta)
        else:
            den = np.maximum(diag - theta, sigma)      # floor sigma: Krylov-like on the low manifold
            t = r / den
        if m == mmax:
            Vn = Y[:, :keep].T @ V[:m]
            AVn = Y[:, :keep].T @ AV[:m]
            V[:keep] = Vn
            AV[:keep] = AVn
            m = keep
        for _ in range(2):
            t = t - V[:m].T @ (V[:m].conj() @ t)
        tn = np.linalg.norm(t)
        if tn < 1e-14:
            t = r - V[:m].T @ (V[:m].conj() @ r)
            tn = np.linalg.norm(t)
            if tn < 1e-14:
                break
        V[m] = t / tn
        AV[m] = mv(V[m])
        m += 1
        nmv += 1
    return theta, u / np.linalg.norm(u), nmv, rn


# ----------------------------------------------------------------------------------------------------------
# Truncation
# ----------------------------------------------------------------------------------------------------------
def _theta_block(th, n0a, Da, n0c, Dc, r):
    """Rows (qa = r, s1 = 0), (qa = 1 - r, s1 = 1); columns (s2 = 0, qc = r), (s2 = 1, qc = 1 - r)."""
    a0, a1 = bsl(n0a, Da, r), bsl(n0a, Da, 1 - r)
    c0, c1 = bsl(n0c, Dc, r), bsl(n0c, Dc, 1 - r)
    top = np.hstack([th[a0, 0, 0, c0], th[a0, 0, 1, c1]])
    bot = np.hstack([th[a1, 1, 0, c0], th[a1, 1, 1, c1]])
    return np.vstack([top, bot])


def split(th, n0a, Da, n0c, Dc, chi, cutoff, direction, alpha=0.0, ts=None, x=None):
    """Truncated split of the two-site tensor.  direction 'right' returns (A left-orthonormal, B carrying the
    weights); 'left' returns (A carrying the weights, B right-orthonormal).  Returns A, B, new n0, discarded
    weight, Schmidt probabilities (from theta alone)."""
    blk = {r: _theta_block(th, n0a, Da, n0c, Dc, r) for r in (0, 1)}
    cand = []          # (weight, r, index)
    dec = {}
    for r in (0, 1):
        Bm = blk[r]
        if Bm.size == 0:
            dec[r] = None
            continue
        if alpha > 0:
            if direction == 'right':
                rho = Bm @ Bm.conj().T
                pert = ts.mixer_rho_left(x, r)
            else:
                rho = Bm.conj().T @ Bm
                pert = ts.mixer_rho_right(x, r)
            tp = np.trace(pert).real
            rho = rho + (alpha * pert / tp if tp > 1e-300 else 0)
            ev, U = np.linalg.eigh((rho + rho.conj().T) / 2)
            ev, U = ev[::-1], U[:, ::-1]
            dec[r] = ('eig', ev, U)
            for i, e in enumerate(ev):
                cand.append((max(e, 0.0), r, i))
        else:
            U, s, Vh = np.linalg.svd(Bm, full_matrices=False)
            dec[r] = ('svd', s, U, Vh)
            for i, e in enumerate(s):
                cand.append((e * e, r, i))
    cand.sort(key=lambda t: -t[0])
    tot = sum(c[0] for c in cand)
    nk = min(chi, len(cand))
    while nk > 1 and cand[nk - 1][0] <= cutoff * cand[0][0]:
        nk -= 1
    kept = {0: [], 1: []}
    for c in cand[:nk]:
        kept[c[1]].append(c[2])
    for r in (0, 1):
        kept[r].sort()
    k0, k1 = len(kept[0]), len(kept[1])
    Dn = k0 + k1
    dt = th.dtype
    A = np.zeros((Da, 2, Dn), dtype=dt)
    B = np.zeros((Dn, 2, Dc), dtype=dt)
    probs = []
    disc = 0.0
    nrm2 = 0.0
    for r in (0, 1):
        if not kept[r]:
            continue
        idx = kept[r]
        rr = slice(0, k0) if r == 0 else slice(k0, Dn)
        na0 = bdim(n0a, Da, r)
        nc0 = bdim(n0c, Dc, r)
        Bm = blk[r]
        if dec[r][0] == 'svd':
            _, s, U, Vh = dec[r]
            Uk, sk, Vk = U[:, idx], s[idx], Vh[idx]
            if direction == 'right':
                L_, R_ = Uk, sk[:, None] * Vk
            else:
                L_, R_ = Uk * sk[None, :], Vk
        else:
            _, ev, U = dec[r]
            if direction == 'right':
                L_ = U[:, idx]
                R_ = L_.conj().T @ Bm
            else:
                R_ = U[:, idx].conj().T          # rows: kept right states
                L_ = Bm @ R_.conj().T
        A[bsl(n0a, Da, r), 0, rr] = L_[:na0]
        A[bsl(n0a, Da, 1 - r), 1, rr] = L_[na0:]
        B[rr, 0, bsl(n0c, Dc, r)] = R_[:, :nc0]
        B[rr, 1, bsl(n0c, Dc, 1 - r)] = R_[:, nc0:]
    # Schmidt weights of the kept (projected) state, and the discarded weight of theta
    nrm2 = np.linalg.norm(A if direction == 'left' else B) ** 2
    tot_th = np.linalg.norm(th) ** 2
    disc = max(0.0, 1.0 - nrm2 / tot_th)
    if direction == 'right':
        B /= np.sqrt(nrm2)
    else:
        A /= np.sqrt(nrm2)
    return A, B, k0, disc


def schmidt_probs_center(C, n0a, Da, n0b, Db):
    """Schmidt probabilities across the right bond of a centre tensor C[a, s, b] (rows (a, s), cols b)."""
    ps = []
    for qb in (0, 1):
        cols = bsl(n0b, Db, qb)
        if bdim(n0b, Db, qb) == 0:
            continue
        Mrow = np.vstack([C[bsl(n0a, Da, qb), 0, cols], C[bsl(n0a, Da, 1 - qb), 1, cols]])
        if Mrow.size:
            ps.append(np.linalg.svd(Mrow, compute_uv=False) ** 2)
    p = np.concatenate(ps)
    return p / p.sum()


# ----------------------------------------------------------------------------------------------------------
# DMRG driver
# ----------------------------------------------------------------------------------------------------------
def sbar(S, L, W):
    S = np.asarray(S)
    l = np.arange(1, L)
    sel = (l >= L // 2 - W) & (l <= L // 2 + W)
    return float(S[sel].mean())


def dmrg(psi, mpo, sweeps, cutoff=1e-14, log=None, ckpt=None, stop=None, energy_shift_per_site=0.0,
         meas_every=None):
    """sweeps: list of dicts {chi, alpha, res_tol, maxit}.  One entry = one right sweep plus one left sweep.
    stop(hist) -> bool may end the run early (checked after each sweep with alpha = 0).
    Returns psi (right-canonical, centre at site 0) and the per-sweep history."""
    L = psi.L
    dt = np.result_type(mpo.dtype, psi.M[0].dtype)
    psi.M = [m.astype(dt) for m in psi.M]
    right_canonicalize(psi)
    W = mpo.W
    Lenv = [None] * (L + 1)
    Renv = [None] * (L + 1)
    w = mpo.w
    Lenv[0] = np.zeros((1, w, 1), dtype=dt); Lenv[0][0, 0, 0] = 1.0
    Renv[L] = np.zeros((1, w, 1), dtype=dt); Renv[L][0, 2, 0] = 1.0
    for j in range(L - 1, 0, -1):
        Renv[j] = env_right(Renv[j + 1], psi.M[j], W[j], psi.n0[j], psi.n0[j + 1], mpo)
    hist = []
    Scut = np.zeros(L - 1)
    Zprof = np.zeros(L)
    t_start = time.time()

    def local(j, sw):
        D = psi.dims
        ts = TwoSite(Lenv[j], W[j], W[j + 1], Renv[j + 2], psi.n0[j], D[j], psi.n0[j + 2], D[j + 2], mpo)
        th = np.tensordot(psi.M[j], psi.M[j + 1], axes=([2], [0]))
        x0 = ts.pack_dense(th)
        if ts.size <= 128:
            H = np.array([ts.matvec(e) for e in np.eye(ts.size, dtype=dt)]).T
            ev, vv = np.linalg.eigh((H + H.conj().T) / 2)
            e, v = ev[0], vv[:, 0]
            r = np.linalg.norm(ts.matvec(v) - e * v)
            nmv = ts.size
        else:
            if not np.any(x0):
                x0 = np.random.default_rng(j).standard_normal(ts.size).astype(dt)
            prec = None
            if sw.get('prec', 'kron4') == 'kron4':
                ts.build_kron4()
                sg = sw.get('sigma', 1e-2)
                prec = lambda rr, thv: ts.precond4(rr, thv, sg)
            e, v, nmv, r = davidson(ts.matvec, x0, ts.diag, res_tol=sw['res_tol'], maxit=sw.get('maxit', 400),
                                    mmax=sw.get('mmax', 32), keep=sw.get('keep', 6), sigma=sw.get('sigma', 1e-3),
                                    prec=prec)
        return ts, e, v, nmv, r

    for isw, sw in enumerate(sweeps):
        t0 = time.time()
        chi, alpha = sw['chi'], sw.get('alpha', 0.0)
        trunc, rmax, nmvs, E = 0.0, 0.0, [], None
        for j in range(L - 1):                                   # right sweep
            ts, E, v, nmv, r = local(j, sw)
            nmvs.append(nmv); rmax = max(rmax, r)
            th = ts.unpack_dense(v)
            A, B, k0, disc = split(th, psi.n0[j], psi.dims[j], psi.n0[j + 2], psi.dims[j + 2], chi, cutoff,
                                   'right', alpha=alpha, ts=ts, x=v)
            trunc = max(trunc, disc)
            psi.M[j], psi.M[j + 1] = A, B
            psi.n0[j + 1] = k0
            Lenv[j + 1] = env_left(Lenv[j], A, W[j], psi.n0[j], psi.n0[j + 1], mpo)
            Renv[j + 2] = None if j + 2 < L else Renv[j + 2]
        for j in range(L - 2, -1, -1):                           # left sweep
            ts, E, v, nmv, r = local(j, sw)
            nmvs.append(nmv); rmax = max(rmax, r)
            th = ts.unpack_dense(v)
            ab = np.abs(th) ** 2
            p1 = ab.sum(axis=(0, 2, 3))
            Zprof[j] = p1[0] - p1[1]
            if j == L - 2:
                p2 = ab.sum(axis=(0, 1, 3))
                Zprof[L - 1] = p2[0] - p2[1]
            A, B, k0, disc = split(th, psi.n0[j], psi.dims[j], psi.n0[j + 2], psi.dims[j + 2], chi, cutoff,
                                   'left', alpha=alpha, ts=ts, x=v)
            trunc = max(trunc, disc)
            psi.M[j], psi.M[j + 1] = A, B
            psi.n0[j + 1] = k0
            Renv[j + 1] = env_right(Renv[j + 2], B, W[j + 1], psi.n0[j + 1], psi.n0[j + 2], mpo)
            if j + 1 < L:
                Lenv[j + 1] = None
            if alpha == 0.0:
                Scut[j] = entropy(schmidt_probs_center(psi.M[j], psi.n0[j], psi.dims[j], psi.n0[j + 1],
                                                       psi.dims[j + 1]))
        rec = dict(sweep=isw, chi=chi, alpha=alpha, E=float(np.real(E)), trunc=trunc, res_max=float(rmax),
                   res_tol=sw['res_tol'], mv_mean=float(np.mean(nmvs)), mv_max=int(np.max(nmvs)),
                   Dmax=int(max(psi.dims)), t=time.time() - t0, t_tot=time.time() - t_start)
        dens = (1 - Zprof) / 2
        rec['n_mid'] = float(np.mean(dens[L // 4: 3 * L // 4]))
        rec['N_tot'] = float(np.sum(dens))
        rec['Zprof'] = Zprof.tolist()
        if alpha == 0.0:
            rec['Sbar1'] = sbar(Scut, L, 1); rec['Sbar6'] = sbar(Scut, L, 6); rec['Sbar12'] = sbar(Scut, L, 12)
            rec['Sprof'] = Scut.tolist()
        if meas_every and alpha == 0.0:
            rec.update(meas_every(psi))
        hist.append(rec)
        if log:
            log(rec)
        if ckpt:
            ckpt(psi, hist)
        if stop is not None and alpha == 0.0 and stop(hist):
            break
    return psi, hist


# ----------------------------------------------------------------------------------------------------------
# Measurements on a right-canonical MPS (centre at site 0)
# ----------------------------------------------------------------------------------------------------------
def _apply_site(E, A, O, Bc=None):
    """E'[b', b] = sum conj(B[a', s', b']) E[a', a] O[s', s] A[a, s, b] (B = A if not given)."""
    if Bc is None:
        Bc = A
    T = np.tensordot(E, A, axes=([1], [0]))                     # a' s b
    if O is not None:
        T = np.tensordot(T, O, axes=([1], [1]))                  # a' b s'
        T = T.transpose(0, 2, 1)
    return np.tensordot(Bc.conj(), T, axes=([0, 1], [0, 1]))     # b' b


def measure(psi, mpo, gph=None, gwin=48, cp_i0=None, cp_rmax=None, h2_list=None, h1_list=None, timing=None):
    """One left-to-right pass.  Returns entropies at all cuts, <Z_j>, bond observables (P_j = XX + YY and
    the bond energies), the boson and fermion one-body matrices on a central window of gwin sites, and the
    connected valley-polarization correlator C_P(r) = <P_i0 P_i0+r> - <P_i0><P_i0+r>.  All operators are in
    the original (unrotated) frame: gph are the site phases of the gauge rotation that was applied to H."""
    t0 = time.time()
    L = psi.L
    psi = psi.copy()
    if gph is None:
        gph = np.ones(L, dtype=complex)
    dt = np.result_type(psi.M[0].dtype, complex if np.any(np.abs(gph - 1) > 1e-14) else float)
    psi.M = [m.astype(dt) for m in psi.M]
    sp = [gph[j] * SP for j in range(L)]                       # U Sp U^dag
    sm = [np.conj(gph[j]) * SM for j in range(L)]
    Zexp = np.zeros(L)
    Pb = np.zeros(L - 1)
    eb = np.zeros(L - 1)
    i_lo = L // 2 - gwin // 2
    win = list(range(i_lo, i_lo + gwin))
    Gb = np.zeros((gwin, gwin), dtype=complex)
    Gf = np.zeros((gwin, gwin), dtype=complex)
    if cp_i0 is None:
        cp_i0 = L // 4
    if cp_rmax is None:
        cp_rmax = L // 2
    CP = None
    for j in range(L):
        C = psi.M[j]                                           # centre tensor (left of it left-canonical)
        rho1 = np.tensordot(C.conj(), C, axes=([0, 2], [0, 2]))  # (s', s)
        Zexp[j] = float(np.real(np.trace(rho1 @ ZZ_)))
        if j < L - 1:
            th = np.tensordot(C, psi.M[j + 1], axes=([2], [0]))  # a s1 s2 c
            rho2 = np.tensordot(th.conj(), th, axes=([0, 3], [0, 3]))  # s1' s2' s1 s2
            rho2 = rho2.reshape(4, 4)
            Pop = 2 * (np.kron(sp[j], sm[j + 1]) + np.kron(sm[j], sp[j + 1]))
            Pb[j] = float(np.real(np.trace(rho2 @ Pop.T)))       # Tr(rho O) with rho[s', s] = <s'|..|s>
            if h2_list is not None:
                h = np.asarray(h2_list[j]) + 0.5 * (np.kron(h1_list[j], I2) + np.kron(I2, h1_list[j + 1]))
                eb[j] = float(np.real(np.trace(rho2 @ np.asarray(h).T)))
        # one-body rows: centre at i in the window
        if j in win:
            ii = j - i_lo
            for lab, G, string in (('b', Gb, I2), ('f', Gf, ZZ_)):
                # i = j: diagonal
                G[ii, ii] = (1 - Zexp[j]) / 2
                T = np.tensordot(C, sm[j], axes=([1], [1]))            # a b s'   (O[s', s] C[a, s, b])
                T = T.transpose(0, 2, 1)                                # a s' b
                E = np.tensordot(C.conj(), T, axes=([0, 1], [0, 1]))    # b'(bra) b(ket)
                for jj in range(j + 1, i_lo + gwin):
                    B = psi.M[jj]
                    # close with sp at jj:  <sm_i ... sp_jj>
                    T = np.tensordot(E, B, axes=([1], [0]))             # b' s c
                    T2 = np.tensordot(T, sp[jj], axes=([1], [1]))       # b' c s'
                    val = np.tensordot(B.conj(), T2, axes=([0, 2, 1], [0, 1, 2]))
                    G[ii, jj - i_lo] = val
                    G[jj - i_lo, ii] = np.conj(val)
                    E = _apply_site(E, B, string)
        # valley-polarization correlator from i0
        if j == cp_i0 and j < L - 3:
            th = np.tensordot(C, psi.M[j + 1], axes=([2], [0]))  # a s1 s2 c
            Pop = (2 * (np.kron(sp[j], sm[j + 1]) + np.kron(sm[j], sp[j + 1]))).reshape(2, 2, 2, 2)
            T = np.tensordot(th, Pop, axes=([1, 2], [2, 3]))     # a c s1' s2'
            T = T.transpose(0, 2, 3, 1)
            E = np.tensordot(th.conj(), T, axes=([0, 1, 2], [0, 1, 2]))    # c' c
            CP = []
            for k in range(j + 2, min(L - 1, j + cp_rmax + 1)):
                B1, B2 = psi.M[k], psi.M[k + 1]
                Pk = (2 * (np.kron(sp[k], sm[k + 1]) + np.kron(sm[k], sp[k + 1]))).reshape(2, 2, 2, 2)
                t2 = np.tensordot(B1, B2, axes=([2], [0]))
                T = np.tensordot(E, t2, axes=([1], [0]))         # c' s1 s2 d
                T = np.tensordot(T, Pk, axes=([1, 2], [2, 3]))   # c' d s1' s2'
                val = np.tensordot(t2.conj(), T, axes=([0, 1, 2, 3], [0, 2, 3, 1]))
                CP.append((k - j, complex(val)))
                E = _apply_site(E, B1, None)
        if j < L - 1:
            _qr_move_right(psi, j)                     # move the centre right (exact, no truncation)
    out = dict(Z=Zexp.tolist(), dens=((1 - Zexp) / 2).tolist(), P_bond=Pb.tolist(), e_bond=eb.tolist(),
               gwin=[i_lo, gwin], Gb_re=Gb.real.tolist(), Gb_im=Gb.imag.tolist(), Gf_re=Gf.real.tolist(),
               Gf_im=Gf.imag.tolist())
    if CP is not None:
        Pm = Pb
        out['CP'] = [(r, float(np.real(v)) - Pm[cp_i0] * Pm[cp_i0 + r]) for r, v in CP]
        out['CP_i0'] = cp_i0
    out['t_meas'] = time.time() - t0
    return out


def entropies(psi):
    """Exact von Neumann entropy at every cut of a right-canonical MPS (centre at site 0)."""
    psi = psi.copy()
    L = psi.L
    S = np.zeros(L - 1)
    for j in range(L - 1):
        _qr_move_right(psi, j)
        # centre now at j+1; Schmidt values across bond j+1 = singular values of the centre as (left bond) x rest
        C = psi.M[j + 1]
        Dl, _, Dr = C.shape
        n0l, n0r = psi.n0[j + 1], psi.n0[j + 2]
        ps = []
        for ql in (0, 1):
            rows = bsl(n0l, Dl, ql)
            if bdim(n0l, Dl, ql) == 0:
                continue
            A = np.hstack([C[rows, 0, bsl(n0r, Dr, ql)], C[rows, 1, bsl(n0r, Dr, 1 - ql)]])
            if A.size:
                ps.append(np.linalg.svd(A, compute_uv=False) ** 2)
        p = np.concatenate(ps)
        S[j] = entropy(p / p.sum())
    return S


def parity_expectation(psi):
    E = np.ones((1, 1))
    for A in psi.M:
        E = _apply_site(E, A, ZZ_)
    return float(np.real(E[0, 0]))


def energy(psi, mpo):
    L = psi.L
    w = mpo.w
    dt = np.result_type(mpo.dtype, psi.M[0].dtype)
    E = np.zeros((1, w, 1), dtype=dt); E[0, 0, 0] = 1.0
    for j in range(L):
        E = env_left(E, psi.M[j].astype(dt), mpo.W[j], psi.n0[j], psi.n0[j + 1], mpo)
    return float(np.real(E[0, 2, 0]))


def variance(psi, mpo):
    """<H^2> - <H>^2 for a normalized MPS, by a double-layer left environment (dense)."""
    L = psi.L
    w = mpo.w
    dt = np.result_type(mpo.dtype, psi.M[0].dtype)
    E = np.zeros((1, w, w, 1), dtype=dt); E[0, 0, 0, 0] = 1.0
    E1 = np.zeros((1, w, 1), dtype=dt); E1[0, 0, 0] = 1.0
    for j in range(L):
        A = psi.M[j].astype(dt)
        W = mpo.W[j]
        T = np.tensordot(E, A, axes=([3], [0]))                       # a' w1 w2 s b
        T = np.tensordot(T, W, axes=([2, 3], [0, 3]))                 # a' w1 b v2 s''
        T = np.tensordot(T, W, axes=([1, 4], [0, 3]))                 # a' b v2 v1 s'
        E = np.tensordot(A.conj(), T, axes=([0, 1], [0, 4]))          # b' b v2 v1
        E = E.transpose(0, 3, 2, 1)
        E1 = env_left(E1, A, W, psi.n0[j], psi.n0[j + 1], mpo)
    h2 = float(np.real(E[0, 2, 2, 0]))
    h1 = float(np.real(E1[0, 2, 0]))
    return h2 - h1 * h1, h1

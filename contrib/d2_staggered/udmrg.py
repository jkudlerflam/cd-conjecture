#!/usr/bin/env python3
"""U(1)-resolved finite two-site DMRG (numpy/scipy only), open or periodic chains.

Every bond index carries an integer charge (number of particles to its left), bond bases are kept
sorted by charge, and the two-site effective Hamiltonian is applied block by block.  The
Hamiltonian conserves the charge exactly, so each run lives in one charge sector N.  Excited states
in a sector are obtained by projecting out previously converged MPS (penalty w |phi><phi|, with
exact overlap environments).

Models (rotated frame of C63 eq. (4), Pauli normalization):
  stag  H' = sum_j [X_j X_j+1 + Y_j Y_j+1 - u Z_j Z_j+1] + b sum_j (-1)^j Z_j, n_j = (1 - Z_j)/2.
        This is unitarily equivalent (alternate-site X rotation, a product unitary) to
        H = sum_j [X X - Y Y + u Z Z + b Z] with the staggered charge Q = sum (-1)^j Z_j = L - 2N.
  tfi   H = -sum X_j X_j+1 - sum Z_j (critical Ising, c = 1/2, v = 2); charge is trivial.
  isxx  local dimension 4: (critical Ising) (x) (XX) on every site, charge = XX particle number.
Written for contrib/d2_staggered (Claude Opus 5.5 subagent, 2026-10-28).
"""
import numpy as np
from scipy.sparse.linalg import eigsh, LinearOperator

X = np.array([[0., 1.], [1., 0.]])
Zp = np.diag([1., -1.])
SP = np.array([[0., 1.], [0., 0.]])   # |0><1| : lowers n (n = 1 for state 1)
SM = SP.T.copy()
I2 = np.eye(2)


# ------------------------------------------------------------------------------------- MPO
def build_mpo(L, d, bonds, onsite, pbc=False):
    """bonds: list of (coef, A, B) for sum_j coef A_j B_j+1 (same for every bond, incl. wrap if pbc).
    onsite: list of L (d,d) arrays.  Left boundary = channel 0, right boundary = channel D-1.
    Returns list of W[i] with shape (Dl, Dr, d_out, d_in)."""
    K = len(bonds)
    nw = K if pbc else 0
    D = 2 + K + nw
    Ws = []
    for i in range(L):
        W = np.zeros((D, D, d, d))
        W[0, 0] = np.eye(d)
        W[D - 1, D - 1] = np.eye(d)
        W[0, D - 1] = onsite[i]
        for k, (c, A, B) in enumerate(bonds):
            W[0, 1 + k] = c * A
            W[1 + k, D - 1] = B
        if pbc:
            for k, (c, A, B) in enumerate(bonds):
                ch = 1 + K + k
                if i == 0:
                    W[0, ch] = B          # bond (L-1, 0): A at L-1, B at 0
                elif i == L - 1:
                    W[ch, D - 1] = c * A
                else:
                    W[ch, ch] = np.eye(d)
        if i == 0:
            W = W[0:1]
        if i == L - 1:
            W = W[:, D - 1:D]
        Ws.append(W)
    return Ws


def model(name, L, u=0.0, b=0.0, pbc=False):
    """Returns (Ws, d, ncharge) with ncharge[s] the particle number of local state s."""
    if name == 'stag':
        bonds = [(2.0, SP, SM), (2.0, SM, SP), (-u, Zp, Zp)]
        onsite = [b * (-1) ** j * Zp for j in range(L)]
        return build_mpo(L, 2, bonds, onsite, pbc), 2, np.array([0, 1])
    if name == 'tfi':
        bonds = [(-1.0, X, X)]
        onsite = [-Zp for j in range(L)]
        return build_mpo(L, 2, bonds, onsite, pbc), 2, np.array([0, 0])
    if name == 'isxx':
        bonds = [(-1.0, np.kron(X, I2), np.kron(X, I2)), (2.0, np.kron(I2, SP), np.kron(I2, SM)),
                 (2.0, np.kron(I2, SM), np.kron(I2, SP))]
        onsite = [-np.kron(Zp, I2) for j in range(L)]
        return build_mpo(L, 4, bonds, onsite, pbc), 4, np.array([0, 1, 0, 1])
    raise KeyError(name)


def channel_charges(Ws, n):
    """delta of each MPO channel (bra minus ket charge), propagated from the left."""
    out = [np.array([0])]
    for W in Ws:
        dl = out[-1]
        dr = np.full(W.shape[1], np.iinfo(np.int64).min)
        for a, b, s, t in zip(*np.nonzero(np.abs(W) > 0)):
            v = dl[a] + n[s] - n[t]
            assert dr[b] in (v, np.iinfo(np.int64).min), 'MPO does not conserve the charge'
            dr[b] = v
        dr[dr == np.iinfo(np.int64).min] = 0
        out.append(dr)
    return out


# ------------------------------------------------------------------------------- utilities
def sectors(q):
    """q sorted -> dict charge -> slice."""
    out = {}
    if len(q) == 0:
        return out
    cut = np.flatnonzero(np.diff(q)) + 1
    st = np.concatenate([[0], cut]); en = np.concatenate([cut, [len(q)]])
    for a, b in zip(st, en):
        out[int(q[a])] = slice(int(a), int(b))
    return out


def block_svd(M, rq, cq, chi, cutoff):
    """M (R, C) block diagonal in rq==cq.  Returns U (R,k), S (k), V (k,C), charges (k) sorted."""
    Us, Ss, Vs, Qs, Ri, Ci = [], [], [], [], [], []
    for q in np.intersect1d(np.unique(rq), np.unique(cq)):
        r = np.flatnonzero(rq == q); c = np.flatnonzero(cq == q)
        blk = M[np.ix_(r, c)]
        if not np.any(blk):
            continue
        try:
            u, s, v = np.linalg.svd(blk, full_matrices=False)
        except np.linalg.LinAlgError:
            import scipy.linalg as sl
            u, s, v = sl.svd(blk, full_matrices=False, lapack_driver='gesvd')
        Us.append(u); Ss.append(s); Vs.append(v); Qs.append(np.full(len(s), q)); Ri.append(r); Ci.append(c)
    allS = np.concatenate(Ss); allQ = np.concatenate(Qs)
    blk_id = np.concatenate([np.full(len(s), i) for i, s in enumerate(Ss)])
    loc = np.concatenate([np.arange(len(s)) for s in Ss])
    order = np.argsort(-allS)
    nrm = np.sum(allS ** 2)
    w = np.cumsum(allS[order] ** 2) / nrm
    keep = min(chi, int(np.sum(allS > 1e-14 * allS[order[0]])))
    # discard tail with weight below cutoff
    k2 = int(np.searchsorted(w, 1 - cutoff)) + 1
    keep = max(1, min(keep, k2))
    sel = order[:keep]
    trunc = float(1 - np.sum(allS[sel] ** 2) / nrm)
    sel = sel[np.lexsort((-allS[sel], allQ[sel]))]
    U = np.zeros((M.shape[0], keep)); V = np.zeros((keep, M.shape[1]))
    for j, idx in enumerate(sel):
        bi, li = blk_id[idx], loc[idx]
        U[Ri[bi], j] = Us[bi][:, li]
        V[j, Ci[bi]] = Vs[bi][li, :]
    S = allS[sel]
    return U, S / np.sqrt(np.sum(S ** 2)), V, allQ[sel], trunc


def lanczos(mv, x0, tol, mmax=30, restarts=6):
    """Lowest eigenpair by restarted Lanczos with full reorthogonalization.
    Stops when the residual norm ||H x - E x|| < tol (absolute)."""
    x = x0 / np.linalg.norm(x0)
    for _ in range(restarts):
        V = [x]; a = []; bb = []
        w = mv(x)
        for m in range(mmax):
            al = V[-1] @ w; a.append(al)
            w = w - al * V[-1] - (bb[-1] * V[-2] if m > 0 else 0)
            Vm = np.array(V)
            w -= Vm.T @ (Vm @ w)
            be = np.linalg.norm(w)
            T = np.diag(a) + np.diag(bb, 1) + np.diag(bb, -1)
            e, y = np.linalg.eigh(T)
            res = abs(be * y[-1, 0])
            if res < tol or be < 1e-14 or m == mmax - 1:
                break
            bb.append(be); V.append(w / be)
            w = mv(V[-1])
        x = np.array(V[:len(y)]).T @ y[:, 0]
        x /= np.linalg.norm(x)
        if res < tol or be < 1e-14:
            break
    return e[0], x


# ------------------------------------------------------------------------------- DMRG
class UDMRG:
    def __init__(self, Ws, d, n, N, seed=1, width=3, mult=2, ortho=(), weight=20.0):
        self.W = Ws; self.L = L = len(Ws); self.d = d; self.n = np.asarray(n); self.N = N
        self.dch = channel_charges(Ws, self.n)
        rng = np.random.default_rng(seed)
        # random initial MPS with charges in a window around the mean filling
        qs = []
        nmax = int(self.n.max())
        for i in range(L + 1):
            lo = max(0, N - nmax * (L - i)); hi = min(nmax * i, N)
            c = int(round(i * N / L))
            win = np.arange(max(lo, c - width), min(hi, c + width) + 1)
            if len(win) == 0:
                win = np.array([min(max(c, lo), hi)])
            m = 1 if i in (0, L) else mult
            qs.append(np.repeat(win, m))
        self.A = []
        for i in range(L):
            T = rng.normal(size=(len(qs[i]), d, len(qs[i + 1])))
            T *= (qs[i][:, None, None] + self.n[None, :, None] == qs[i + 1][None, None, :])
            self.A.append(T)
        self.q = qs
        self.ortho = list(ortho)    # previously converged states (UDMRG objects, canonical)
        self.weight = weight
        self._right_canonicalize()
        self.LE = [None] * (L + 1); self.RE = [None] * (L + 1)
        self.LE[0] = np.ones((1, 1, 1)); self.RE[L] = np.ones((1, 1, 1))
        self.LO = [[None] * (L + 1) for _ in self.ortho]; self.RO = [[None] * (L + 1) for _ in self.ortho]
        for k in range(len(self.ortho)):
            self.LO[k][0] = np.ones((1, 1)); self.RO[k][L] = np.ones((1, 1))
        for i in range(L - 1, 0, -1):
            self._update_R(i)
        self.S = [None] * (L + 1); self.Sq = [None] * (L + 1)
        self.trunc = 0.0

    def _right_canonicalize(self):
        L = self.L
        for i in range(L - 1, 0, -1):
            T = self.A[i]
            chiL, d, chiR = T.shape
            M = T.reshape(chiL, d * chiR)
            cq = (self.q[i + 1][None, :] - self.n[:, None]).ravel()
            U, S, V, nq, _ = block_svd(M, self.q[i], cq, 10 ** 9, 0.0)
            self.A[i] = V.reshape(-1, d, chiR)
            # U rows are the old bond basis; fold U S into the left tensor
            self.A[i - 1] = np.tensordot(self.A[i - 1], U * S[None, :], ([2], [0]))
            self.q[i] = nq
        nrm = np.linalg.norm(self.A[0]); self.A[0] /= nrm

    # environments --------------------------------------------------------------------
    def _update_L(self, i):
        A = self.A[i]
        T = np.tensordot(self.LE[i], A, ([2], [0]))            # a,w,s',b'
        T = np.tensordot(T, self.W[i], ([1, 2], [0, 3]))       # a,b',v,s
        T = np.tensordot(A, T, ([0, 1], [0, 3]))               # b,b',v
        self.LE[i + 1] = T.transpose(0, 2, 1)
        for k, ph in enumerate(self.ortho):
            T = np.tensordot(self.LO[k][i], A, ([1], [0]))      # a',s,b
            self.LO[k][i + 1] = np.tensordot(ph.A[i], T, ([0, 1], [0, 1]))

    def _update_R(self, i):
        B = self.A[i]
        T = np.tensordot(B, self.RE[i + 1], ([2], [2]))       # a',s',b,y
        T = np.tensordot(T, self.W[i], ([1, 3], [3, 1]))       # a',b,w,s
        T = np.tensordot(B, T, ([1, 2], [3, 1]))               # a,a',w
        self.RE[i] = T.transpose(0, 2, 1)
        for k, ph in enumerate(self.ortho):
            T = np.tensordot(B, self.RO[k][i + 1], ([2], [1]))  # a,s,b'
            self.RO[k][i] = np.tensordot(T, ph.A[i], ([1, 2], [1, 2])).T

    # two-site effective Hamiltonian ----------------------------------------------------
    def _heff(self, i):
        Le, Re = self.LE[i], self.RE[i + 2]
        W1, W2 = self.W[i], self.W[i + 1]
        W12 = np.tensordot(W1, W2, ([1], [0]))                 # w,s,s',y,t,t'
        W12 = W12.transpose(0, 3, 1, 4, 2, 5)                  # w,y,s,t,s',t'
        nz = list(zip(*np.nonzero(np.abs(W12) > 1e-15)))
        qL, qR = self.q[i], self.q[i + 2]
        sL, sR = sectors(qL), sectors(qR)
        dl, dr = self.dch[i], self.dch[i + 2]
        n = self.n; d = self.d
        mask = (qL[:, None, None, None] + n[None, :, None, None] + n[None, None, :, None]
                == qR[None, None, None, :])
        idx = np.flatnonzero(mask.ravel())
        shape = (len(qL), d, d, len(qR))
        # group W12 entries by (w, s', t')
        terms = {}
        for (w, y, s, t, sp, tp) in nz:
            terms.setdefault((w, sp, tp), []).append((y, s, t, W12[w, y, s, t, sp, tp]))
        in_blocks = []
        for qa, sla in sL.items():
            for sp in range(d):
                for tp in range(d):
                    qb = qa + n[sp] + n[tp]
                    if qb in sR:
                        in_blocks.append((qa, sla, sp, tp, qb, sR[qb]))

        # precomputed contraction plan
        keys = {}; outs = []
        plan = []
        for (qa_, sla_, sp, tp, qb_, slb_) in in_blocks:
            for w in range(Le.shape[1]):
                lst = terms.get((w, sp, tp))
                if not lst:
                    continue
                qa = qa_ + dl[w]
                if qa not in sL:
                    continue
                ent = []
                for (y, s, t, cf) in lst:
                    qb = qa + n[s] + n[t]
                    qbp = qb - dr[y]
                    if qb not in sR or qbp not in sR:
                        continue
                    key = (y, qa, s, t)
                    if key not in keys:
                        keys[key] = len(outs)
                        outs.append((sL[qa], s, t, sR[qb], np.ascontiguousarray(Re[sR[qb], y, sR[qbp]].T)))
                    ent.append((keys[key], cf))
                if ent:
                    plan.append((sla_, sp, tp, slb_, np.ascontiguousarray(Le[sL[qa], w, sla_]), ent))

        def apply(theta):
            out = np.zeros(shape)
            Z = [None] * len(outs)
            for (sla_, sp, tp, slb_, Lb, ent) in plan:
                Y = Lb @ theta[sla_, sp, tp, slb_]
                for (k, cf) in ent:
                    if Z[k] is None:
                        Z[k] = cf * Y
                    else:
                        Z[k] += cf * Y
            for k, (sa, s, t, sb, Rb) in enumerate(outs):
                if Z[k] is not None:
                    out[sa, s, t, sb] += Z[k] @ Rb
            return out
        return apply, idx, shape

    def _projectors(self, i):
        ps = []
        for k, ph in enumerate(self.ortho):
            T = np.tensordot(self.LO[k][i], ph.A[i], ([0], [0]))       # a,s,c
            T = np.tensordot(T, ph.A[i + 1], ([2], [0]))                # a,s,t,b'
            T = np.tensordot(T, self.RO[k][i + 2], ([3], [0]))          # a,s,t,b
            ps.append(T)
        return ps

    def _optimize(self, i, tol):
        apply, idx, shape = self._heff(i)
        ps = [p.ravel()[idx] for p in self._projectors(i)]
        theta0 = np.tensordot(self.A[i], self.A[i + 1], ([2], [0])).ravel()[idx]
        dim = len(idx)
        w = self.weight

        def mv(x):
            full = np.zeros(int(np.prod(shape))); full[idx] = x
            y = apply(full.reshape(shape)).ravel()[idx]
            for p in ps:
                y += w * p * (p @ x)
            return y
        if dim <= 300:
            H = np.column_stack([mv(e) for e in np.eye(dim)])
            H = (H + H.T) / 2
            e, v = np.linalg.eigh(H); E, x = e[0], v[:, 0]
        else:
            if not np.any(theta0):
                theta0 = np.random.default_rng(0).normal(size=dim)
            E, x = lanczos(mv, theta0, tol)
        ov = sum((p @ x) ** 2 for p in ps)
        Eh = E - w * ov
        full = np.zeros(int(np.prod(shape))); full[idx] = x
        return Eh, full.reshape(shape), float(ov)

    def _split(self, i, theta, chi, cutoff, right, noise=0.0):
        chiL, d, _, chiR = theta.shape
        if noise > 0:
            m = (self.q[i][:, None, None, None] + self.n[None, :, None, None] + self.n[None, None, :, None]
                 == self.q[i + 2][None, None, None, :])
            theta = theta + noise * np.random.default_rng().normal(size=theta.shape) * m
            theta /= np.linalg.norm(theta)
        M = theta.reshape(chiL * d, d * chiR)
        rq = (self.q[i][:, None] + self.n[None, :]).ravel()
        cq = (self.q[i + 2][None, :] - self.n[:, None]).ravel()
        U, S, V, nq, tr = block_svd(M, rq, cq, chi, cutoff)
        self.q[i + 1] = nq
        self.S[i + 1] = S; self.Sq[i + 1] = nq
        self.trunc = max(self.trunc, tr)
        if right:
            self.A[i] = U.reshape(chiL, d, -1)
            self.A[i + 1] = (S[:, None] * V).reshape(-1, d, chiR)
        else:
            self.A[i] = (U * S[None, :]).reshape(chiL, d, -1)
            self.A[i + 1] = V.reshape(-1, d, chiR)

    def sweep(self, chi, cutoff=1e-12, tol=1e-6, noise=0.0):
        """One right sweep then one left sweep; returns energy and max overlap^2 with ortho states."""
        L = self.L; self.trunc = 0.0; ovmax = 0.0
        for i in range(0, L - 1):
            E, th, ov = self._optimize(i, tol)
            ovmax = max(ovmax, ov)
            self._split(i, th, chi, cutoff, right=True, noise=noise)
            self._update_L(i)
        for i in range(L - 2, -1, -1):
            E, th, ov = self._optimize(i, tol)
            ovmax = max(ovmax, ov)
            self._split(i, th, chi, cutoff, right=False, noise=noise)
            self._update_R(i + 1)
        self.E = E
        return E, ovmax

    def measure(self):
        """Exact Schmidt data on every bond from a left-to-right canonical pass (no optimization).
        Leaves the MPS left canonical with the centre on the last site."""
        L = self.L
        out_S, out_F, out_spec = [0.0], [0.0], [None]
        for i in range(L - 1):
            T = self.A[i]; chiL, d, chiR = T.shape
            M = T.reshape(chiL * d, chiR)
            rq = (self.q[i][:, None] + self.n[None, :]).ravel()
            U, S, V, nq, _ = block_svd(M, rq, self.q[i + 1], 10 ** 9, 0.0)
            self.A[i] = U.reshape(chiL, d, -1)
            self.A[i + 1] = np.tensordot(S[:, None] * V, self.A[i + 1], ([1], [0]))
            self.q[i + 1] = nq
            p = S ** 2; p = p[p > 1e-16]; qq = nq[S ** 2 > 1e-16]
            out_S.append(float(-np.sum(p * np.log(p))))
            out_F.append(float(np.sum(p * qq ** 2) - np.sum(p * qq) ** 2))
            spec = {}
            for qv, pv in zip(qq, p):
                spec.setdefault(int(qv), []).append(float(pv))
            out_spec.append({k: v[:8] for k, v in spec.items()})
        out_S.append(0.0); out_F.append(0.0)
        # restore right-canonical form with centre at site 0 for further sweeps
        for i in range(L - 1, 0, -1):
            T = self.A[i]; chiL, d, chiR = T.shape
            M = T.reshape(chiL, d * chiR)
            cq = (self.q[i + 1][None, :] - self.n[:, None]).ravel()
            U, S, V, nq, _ = block_svd(M, self.q[i], cq, 10 ** 9, 0.0)
            self.A[i] = V.reshape(-1, d, chiR)
            self.A[i - 1] = np.tensordot(self.A[i - 1], U * S[None, :], ([2], [0]))
            self.q[i] = nq
        for i in range(L - 1, 0, -1):
            self._update_R(i)
        return out_S, out_F, out_spec

    def energy(self):
        """<psi|H|psi> by a full left contraction."""
        Le = np.ones((1, 1, 1))
        for i in range(self.L):
            A = self.A[i]
            T = np.tensordot(Le, A, ([2], [0])); T = np.tensordot(T, self.W[i], ([1, 2], [0, 3]))
            Le = np.tensordot(A, T, ([0, 1], [0, 3])).transpose(0, 2, 1)
        nrm = sum(np.sum(a * a) for a in [self.A[0]]) if False else None
        return float(Le.ravel()[0])

    def overlap(self, other):
        E = np.ones((1, 1))
        for i in range(self.L):
            T = np.tensordot(E, self.A[i], ([1], [0]))
            E = np.tensordot(other.A[i], T, ([0, 1], [0, 1]))
        return float(E.ravel()[0])

    def norm2(self):
        return self.overlap(self)

    def bond_chis(self):
        return [a.shape[2] for a in self.A[:-1]]


def run_state(Ws, d, n, N, chis, sweeps=3, maxsweeps=12, etol=1e-9, ortho=(), seed=1, log=None, noise0=1e-4):
    """Grow chi through chis; at each chi sweep until |dE| < etol*L (min `sweeps` sweeps)."""
    dm = UDMRG(Ws, d, n, N, seed=seed, ortho=ortho)
    L = len(Ws); hist = []
    for ci, chi in enumerate(chis):
        Eold = None
        for sw in range(maxsweeps):
            nz = noise0 if (ci == 0 and sw < 2) else 0.0
            E, ov = dm.sweep(chi, noise=nz)
            if log:
                log(f'  N={N} chi={chi} sw={sw} E={E:.12f} trunc={dm.trunc:.2e} ov={ov:.1e} maxchi={max(dm.bond_chis())}')
            if sw + 1 >= sweeps and Eold is not None and abs(E - Eold) < etol * L and nz == 0:
                break
            Eold = E
        hist.append(dict(chi=chi, E=float(E), trunc=float(dm.trunc), sweeps=sw + 1, ov=float(ov),
                         maxchi=int(max(dm.bond_chis()))))
    return dm, hist

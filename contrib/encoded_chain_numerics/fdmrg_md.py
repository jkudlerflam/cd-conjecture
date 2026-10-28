#!/usr/bin/env python3
"""Two-site open-chain DMRG with site-dependent dimensions and MPO bonds, plus low excited
states by orthogonal penalty (H + w sum_k |psi_k><psi_k|).  Adapted from
../d2_followup_220_1282/fdmrg.py (same conventions, contractions imported from vumps.py).

usage:
  fdmrg_md.py encoded --q 2 --qp 1 --n 1 --g 0.5 --lam 0.1 --L 32 --chis 32 64 128 [--nstates 2]
  fdmrg_md.py logical --q 2 --qp 1 --n 1 --L 64 --chis ...   (first-order chain, t = 1)
  fdmrg_md.py ising   --L 64 --chis ...
Saves the energies, all bond entropies (cell-boundary bonds are every third bond for the encoded
chain) and cell-resolved bond energies.

Written by Claude (Opus 5.5) for contrib/encoded_chain_numerics, 2026-10-27.
"""
import os, sys, json, time, argparse
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '4')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import eigsh, LinearOperator
from vumps import lt, rt, hac2      # noqa: E402
import cells as C                    # noqa: E402


def ov_l(X, phi, psi):
    t = np.tensordot(X, psi, ([1], [0]))                  # a', s, b
    return np.tensordot(phi.conj(), t, ([0, 1], [0, 1]))  # b', b


def ov_r(X, phi, psi):
    t = np.tensordot(psi, X, ([2], [1]))                  # a, s, b'
    return np.tensordot(phi.conj(), t, ([1, 2], [1, 2]))    # (a', a): bra index first


class FDMRG:
    def __init__(self, Ws, dims, seed=1, chi0=8, ortho=(), weight=10.0, bmult=None):
        self.L = L = len(dims)
        self.bmult = bmult if bmult is not None else [1] * (L - 1)   # chi cap of bond i = bmult[i] * chi
        self.Ws, self.dims = Ws, dims
        rng = np.random.default_rng(seed)
        bd = [1]
        for i in range(1, L):
            left = np.prod(np.array(dims[:i], float)); right = np.prod(np.array(dims[i:], float))
            bd.append(int(min(chi0, left, right)))
        bd.append(1)
        M = [rng.normal(size=(bd[i], dims[i], bd[i + 1])) + 1j * rng.normal(size=(bd[i], dims[i], bd[i + 1]))
             for i in range(L)]
        for i in range(L - 1, 0, -1):
            a, s, b = M[i].shape
            qm, r = np.linalg.qr(M[i].reshape(a, s * b).T)
            M[i] = qm.T.reshape(-1, s, b)
            M[i - 1] = np.tensordot(M[i - 1], r.T, ([2], [0]))
        M[0] /= np.linalg.norm(M[0])
        self.M = M
        D0 = Ws[0].shape[0]
        DL = Ws[-1].shape[1]
        self.Le = [None] * (L + 1)
        self.Re = [None] * (L + 1)
        self.Le[0] = np.zeros((1, D0, 1), complex); self.Le[0][0, 0, 0] = 1
        self.Re[L] = np.zeros((1, DL, 1), complex); self.Re[L][0, DL - 1, 0] = 1
        for i in range(L - 1, 0, -1):
            self.Re[i] = rt(self.Re[i + 1], M[i], Ws[i])
        self.ortho = [list(p) for p in ortho]
        self.w = weight
        self.LO = [[None] * (L + 1) for _ in self.ortho]
        self.RO = [[None] * (L + 1) for _ in self.ortho]
        for k, phi in enumerate(self.ortho):
            self.LO[k][0] = np.ones((1, 1), complex)
            self.RO[k][L] = np.ones((1, 1), complex)
            for i in range(L - 1, 0, -1):
                self.RO[k][i] = ov_r(self.RO[k][i + 1], phi[i], M[i])
        self.S = np.zeros(L - 1)
        self.E = None

    def _proj(self, i):
        ps = []
        for k, phi in enumerate(self.ortho):
            pp = np.tensordot(phi[i], phi[i + 1], ([2], [0]))         # a' s1 s2 b'
            t = np.tensordot(self.LO[k][i], pp.conj(), ([0], [0]))     # a s1 s2 b'
            t = np.tensordot(t, self.RO[k][i + 2], ([3], [0]))         # a s1 s2 b
            ps.append(t.conj())
        return ps

    def sweep(self, chi, tol=1e-10, cutoff=1e-14):
        L, Ws, M = self.L, self.Ws, self.M
        maxtr = 0.0
        E = None
        for direction in (1, -1):
            rng = range(L - 1) if direction == 1 else range(L - 2, -1, -1)
            for i in rng:
                th = np.tensordot(M[i], M[i + 1], ([2], [0]))
                shp = th.shape
                ps = self._proj(i)
                Le, Re, W1, W2 = self.Le[i], self.Re[i + 2], Ws[i], Ws[i + 1]
                w = self.w

                def f(v):
                    x = v.reshape(shp)
                    y = hac2(Le, W1, W2, Re, x)
                    for p in ps:
                        y = y + w * p * np.vdot(p, x)
                    return y.ravel()
                nn = th.size
                if nn <= 300:
                    H = np.zeros((nn, nn), complex)
                    for k in range(nn):
                        e = np.zeros(nn, complex); e[k] = 1
                        H[:, k] = f(e)
                    ev, Uv = np.linalg.eigh((H + H.conj().T) / 2)
                    E, th = ev[0], Uv[:, 0].reshape(shp)
                else:
                    op = LinearOperator((nn, nn), matvec=f, dtype=complex)
                    v0 = th.ravel()
                    ev, Uv = eigsh(op, k=1, which='SA', v0=v0, tol=tol, ncv=min(nn - 1, 24), maxiter=10000)
                    E, th = ev[0], Uv[:, 0].reshape(shp)
                a, s1, s2, b = shp
                u, s, vh = linalg.svd(th.reshape(a * s1, s2 * b), full_matrices=False, lapack_driver='gesdd')
                keep = min(int(chi * self.bmult[i]), max(1, int((s > cutoff * s[0]).sum())))
                maxtr = max(maxtr, float((s[keep:] ** 2).sum()))
                s = s[:keep] / np.linalg.norm(s[:keep])
                p = s ** 2; p = p[p > 1e-30]
                self.S[i] = float(-(p * np.log(p)).sum())
                if direction == 1:
                    M[i] = u[:, :keep].reshape(a, s1, keep)
                    M[i + 1] = (s[:, None] * vh[:keep]).reshape(keep, s2, b)
                    self.Le[i + 1] = lt(self.Le[i], M[i], Ws[i])
                    for k, phi in enumerate(self.ortho):
                        self.LO[k][i + 1] = ov_l(self.LO[k][i], phi[i], M[i])
                else:
                    M[i + 1] = vh[:keep].reshape(keep, s2, b)
                    M[i] = (u[:, :keep] * s[None, :]).reshape(a, s1, keep)
                    self.Re[i + 1] = rt(self.Re[i + 2], M[i + 1], Ws[i + 1])
                    for k, phi in enumerate(self.ortho):
                        self.RO[k][i + 1] = ov_r(self.RO[k][i + 2], phi[i + 1], M[i + 1])
        self.E = float(np.real(E))
        return self.E, maxtr

    def energy(self):
        """<H> without the penalty (state is right-canonical with centre on site 0 after a sweep)."""
        X = self.Le[0]
        for i in range(self.L):
            X = lt(X, self.M[i], self.Ws[i])
        return float(np.real(X[0, -1, 0]))

    def overlaps(self):
        out = []
        for phi in self.ortho:
            X = np.ones((1, 1), complex)
            for i in range(self.L):
                X = ov_l(X, phi[i], self.M[i])
            out.append(abs(X[0, 0]))
        return out

    def two_site(self, i, O):
        """<O> on sites (i, i+1); requires a left-canonical prefix -> compute via norm envs."""
        X = np.ones((1, 1), complex)
        for j in range(i):
            X = ov_l(X, self.M[j], self.M[j])
        th = np.tensordot(self.M[i], self.M[i + 1], ([2], [0]))
        Y = np.ones((1, 1), complex)
        for j in range(self.L - 1, i + 1, -1):
            Y = ov_r(Y, self.M[j], self.M[j])
        a, d1, d2, b = th.shape
        t = np.tensordot(X, th, ([1], [0]))
        t = np.tensordot(t, Y, ([3], [1]))                     # (a', s1, s2, b')
        Ot = np.einsum('pq,aqb->apb', O, t.reshape(a, d1 * d2, b))
        return float(np.real(np.vdot(th.reshape(a, d1 * d2, b), Ot)))


def build(a):
    if a.model == 'ising':
        d, t = C.ising_terms()
        W = C.mpo_from_terms(d, t)
        return [2] * a.L, [W[0]] * a.L, 1, None
    c = C.cell(a.q, a.qp, a.n, 1.0, a.g)
    if a.model == 'logical':
        dims, terms = C.logical_chain_terms(c, 1.0, 1)
        W = C.mpo_from_terms(dims, terms)
        return dims * a.L, W * a.L, 1, c
    dims, W = C.chain_mpo(c, a.lam, 1)
    return dims * a.L, W * a.L, 3, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('model', choices=['ising', 'logical', 'encoded'])
    ap.add_argument('--q', type=int, default=2)
    ap.add_argument('--qp', type=int, default=1)
    ap.add_argument('--n', type=int, default=1)
    ap.add_argument('--g', type=float, default=0.5)
    ap.add_argument('--lam', type=float, default=0.1)
    ap.add_argument('--L', type=int, required=True, help='number of cells (sites for ising)')
    ap.add_argument('--chis', type=int, nargs='+', default=[16, 32, 64])
    ap.add_argument('--sweeps', type=int, default=4)
    ap.add_argument('--maxsweeps', type=int, default=20)
    ap.add_argument('--etol', type=float, default=1e-9)
    ap.add_argument('--nstates', type=int, default=1)
    ap.add_argument('--weight', type=float, default=5.0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', required=True)
    ap.add_argument('--sparse', action='store_true', help='sparse MPO contractions (use for q = 3)')
    a = ap.parse_args()
    dims, Ws, sub, c = build(a)
    if a.sparse:
        import sparse_mpo as SPM
        global lt, rt, hac2
        lt, rt, hac2 = SPM.lt, SPM.rt, SPM.hac2
        cache = {}
        Ws = [cache.setdefault(id(w), SPM.from_dense(w)) for w in Ws]
    data = dict(args=vars(a), sub=sub, states=[])
    t0 = time.time()
    found = []
    for st in range(a.nstates):
        bm = [1 if (i + 1) % sub == 0 else (a.q if sub > 1 else 1) for i in range(len(dims) - 1)]
        dm = FDMRG(Ws, dims, seed=a.seed + st, ortho=found, weight=a.weight, bmult=bm)
        recs = []
        for chi in a.chis:
            Eprev = None
            for sw in range(a.maxsweeps):
                E, tr = dm.sweep(chi)
                if sw + 1 >= a.sweeps and Eprev is not None and abs(E - Eprev) < a.etol * a.L:
                    break
                Eprev = E
            Eh = dm.energy()
            rec = dict(chi=chi, E=Eh, Epen=E, sweeps=sw + 1, trunc=tr, S=dm.S.tolist(), ov=dm.overlaps(),
                       time_s=time.time() - t0)
            if a.model == 'encoded' and chi == a.chis[-1]:
                bond = sum(np.kron(f, f) for f in c['F'])
                rec['e_bond'] = [a.lam * dm.two_site(3 * j + 2, bond) for j in range(a.L - 1)]
            recs.append(rec)
            Sb = dm.S[sub - 1::sub] if sub > 1 else dm.S
            print(f"state {st} L={a.L} chi={chi} E={Eh:.12f} sweeps={sw+1} trunc={tr:.1e} "
                  f"Smid={Sb[len(Sb)//2]:.6f} ov={np.round(rec['ov'], 8).tolist()} t={time.time()-t0:.0f}s", flush=True)
            data['states'] = [s for s in data['states'] if s['state'] != st] + [dict(state=st, recs=recs)]
            json.dump(data, open(a.out + '.tmp', 'w'))
            os.replace(a.out + '.tmp', a.out)
        found.append([m.copy() for m in dm.M])


if __name__ == '__main__':
    main()

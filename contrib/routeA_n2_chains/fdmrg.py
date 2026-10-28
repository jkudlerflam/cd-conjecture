#!/usr/bin/env python3
"""Two-site finite DMRG on an open chain with a general finite-range MPO, plus bulk measurements.

Used here in two ways.
 (a) Bulk finite-entanglement scaling: a long chain (L = 400 to 800) at small chi, where the
     chi-induced correlation length xi_chi is much shorter than the distance to the edges, so the
     central region is the infinite-system finite-chi state.  xi_chi is read from the transfer
     matrix of a window of m central sites, -m / ln|lambda_1|, and S from the central bonds.
 (b) Ordinary open-chain entropy profiles at L = 96 to 384 with large chi.

usage: fdmrg.py MODEL --L 400 --chis 8 12 16 24 32 48 --sweeps 6 [--seed 1] [--win 64 128]
                [--out data/fdmrg_MODEL_L400.json]
Conventions as in vumps.py (MPO channel 0 = start, D-1 = done; A[k] has shape (chi_l, d, chi_r)).
Written for contrib/d2_followup_220_1282 (Claude Opus 5.5, 2026-10-27).
"""
import os, sys
for _k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_k, "4")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import argparse, json, time
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import eigsh, eigs, LinearOperator
from vumps import lt, rt, hac2, mpo_from_term, PAULI, tr_id


def two_site_min(Lenv, W1, W2, Renv, th, tol):
    shp = th.shape
    n = th.size
    f = lambda v: hac2(Lenv, W1, W2, Renv, v.reshape(shp)).ravel()
    if n <= 256:
        H = np.zeros((n, n), complex)
        for k in range(n):
            e = np.zeros(n, complex); e[k] = 1
            H[:, k] = f(e)
        w, U = np.linalg.eigh((H + H.conj().T) / 2)
        return w[0], U[:, 0].reshape(shp)
    op = LinearOperator((n, n), matvec=f, dtype=complex)
    w, U = eigsh(op, k=1, which='SA', v0=th.ravel(), tol=tol, ncv=min(n - 1, 24), maxiter=10000)
    return w[0], U[:, 0].reshape(shp)


class FDMRG:
    def __init__(self, W, L, d=2, seed=1, chi0=8, edge=None):
        self.L, self.d = L, d
        self.Ws = [W] * L
        if edge is not None:                       # extra one-site operator on both end sites
            D0 = W.shape[0]
            for j in (0, L - 1):
                Wj = W.copy(); Wj[0, D0 - 1] = Wj[0, D0 - 1] + edge
                self.Ws[j] = Wj
        D = W.shape[0]
        rng = np.random.default_rng(seed)
        dims = [1] + [min(chi0, d ** i, d ** (L - i)) for i in range(1, L)] + [1]
        M = [rng.normal(size=(dims[i], d, dims[i + 1])) + 1j * rng.normal(size=(dims[i], d, dims[i + 1])) for i in range(L)]
        for i in range(L - 1, 0, -1):            # right-canonicalize
            a, s, b = M[i].shape
            q, r = np.linalg.qr(M[i].reshape(a, s * b).T)
            M[i] = q.T.reshape(-1, s, b)
            M[i - 1] = np.tensordot(M[i - 1], r.T, ([2], [0]))
        M[0] /= np.linalg.norm(M[0])
        self.M = M
        self.Le = [None] * (L + 1)
        self.Re = [None] * (L + 1)
        self.Le[0] = np.zeros((1, D, 1), complex); self.Le[0][0, 0, 0] = 1
        self.Re[L] = np.zeros((1, D, 1), complex); self.Re[L][0, D - 1, 0] = 1
        for i in range(L - 1, 0, -1):
            self.Re[i] = rt(self.Re[i + 1], M[i], self.Ws[i])
        self.S = np.zeros(L - 1)
        self.E = None

    def sweep(self, chi, tol=1e-10, cutoff=1e-14):
        L, Ws, M = self.L, self.Ws, self.M
        maxtr = 0.0
        for direction in (1, -1):
            rng = range(L - 1) if direction == 1 else range(L - 2, -1, -1)
            for i in rng:
                th = np.tensordot(M[i], M[i + 1], ([2], [0]))
                E, th = two_site_min(self.Le[i], Ws[i], Ws[i + 1], self.Re[i + 2], th, tol)
                a, s1, s2, b = th.shape
                u, s, vh = linalg.svd(th.reshape(a * s1, s2 * b), full_matrices=False, lapack_driver='gesdd')
                keep = min(chi, max(1, int((s > cutoff * s[0]).sum())))
                maxtr = max(maxtr, float((s[keep:] ** 2).sum()))
                s = s[:keep] / np.linalg.norm(s[:keep])
                p = s ** 2; p = p[p > 1e-30]
                self.S[i] = float(-(p * np.log(p)).sum())
                if direction == 1:
                    M[i] = u[:, :keep].reshape(a, s1, keep)
                    M[i + 1] = (s[:, None] * vh[:keep]).reshape(keep, s2, b)
                    self.Le[i + 1] = lt(self.Le[i], M[i], Ws[i])
                else:
                    M[i + 1] = vh[:keep].reshape(keep, s2, b)
                    M[i] = (u[:, :keep] * s[None, :]).reshape(a, s1, keep)
                    self.Re[i + 1] = rt(self.Re[i + 2], M[i + 1], Ws[i + 1])
        self.E = float(E)
        return self.E, maxtr

    # ---- measurements (state is right-canonical with the centre on site 0 after a sweep)
    def local_ops(self, ops):
        """<O_j> for every site j and every operator in ops (list of d x d)."""
        M = self.M
        out = np.zeros((len(ops), self.L), complex)
        X = np.eye(1, dtype=complex)
        for j in range(self.L):
            A = M[j]
            t = np.tensordot(X, A, ([1], [0]))                     # a', s, b
            for k, O in enumerate(ops):
                t2 = np.tensordot(t, O, ([1], [1]))                  # a', b, s'
                out[k, j] = np.tensordot(A.conj(), t2, ([0, 1, 2], [0, 2, 1]))
            X = np.tensordot(A.conj(), t, ([0, 1], [0, 1]))
        return out.real

    def left_norm_envs(self):
        X = [np.eye(1, dtype=complex)]
        for j in range(self.L):
            A = self.M[j]
            X.append(np.tensordot(A.conj(), np.tensordot(X[-1], A, ([1], [0])), ([0, 1], [0, 1])))
        return X

    def bond_energies(self, h2):
        """<h_{j,j+1}> for a two-site operator h2 (4 x 4)."""
        X = self.left_norm_envs()
        H = h2.reshape(2, 2, 2, 2)
        out = []
        for j in range(self.L - 1):
            th = np.tensordot(self.M[j], self.M[j + 1], ([2], [0]))   # a s1 s2 b (right part is canonical)
            t = np.tensordot(X[j], th, ([1], [0]))
            Ht = np.einsum('stuv,auvb->astb', H, t)
            out.append(float(np.tensordot(th.conj(), Ht, ([0, 1, 2, 3], [0, 1, 2, 3])).real))
        return np.array(out)

    def correlations(self, i0, rmax, ops):
        """C_ab(r) = <O_a(i0) O_b(i0+r)> - <O_a><O_b>, r = 1..rmax."""
        X = self.left_norm_envs()
        loc = self.local_ops(ops)
        na = len(ops)
        C = np.zeros((na, na, rmax), complex)
        A = self.M[i0]
        for a, Oa in enumerate(ops):
            t = np.tensordot(X[i0], A, ([1], [0]))
            t = np.tensordot(t, Oa, ([1], [1])).transpose(0, 2, 1)       # a', s', b
            Y = np.tensordot(A.conj(), t, ([0, 1], [0, 1]))              # b', b
            for r in range(1, rmax + 1):
                B = self.M[i0 + r]
                tb = np.tensordot(Y, B, ([1], [0]))                       # b', s, c
                for b, Ob in enumerate(ops):
                    t2 = np.tensordot(tb, Ob, ([1], [1]))                  # b', c, s'
                    C[a, b, r - 1] = np.tensordot(B.conj(), t2, ([0, 1, 2], [0, 2, 1]))
                Y = np.tensordot(B.conj(), tb, ([0, 1], [0, 1]))
            for b in range(na):
                C[a, b] -= loc[a, i0] * loc[b, i0 + np.arange(1, rmax + 1)]
        return C

    def window_spectrum(self, j0, m, nev=6):
        """Eigenvalues of prod_{j0 <= j < j0+m} E_{B_j} (right-canonical tensors; E(x) = sum_s B x B^dag)."""
        Bs = self.M[j0:j0 + m]
        chi = Bs[0].shape[0]
        assert Bs[-1].shape[2] == chi
        N = chi * chi

        def mv(v):
            X = v.reshape(chi, chi)
            for B in Bs[::-1]:
                X = tr_id(X, B)
            return X.ravel()
        if N <= 400:
            Mx = np.zeros((N, N), complex)
            for k in range(N):
                e = np.zeros(N, complex); e[k] = 1
                Mx[:, k] = mv(e)
            w = linalg.eigvals(Mx)
        else:
            w = eigs(LinearOperator((N, N), matvec=mv, dtype=complex), k=nev, which='LM', tol=1e-10,
                     return_eigenvectors=False, maxiter=5000, ncv=max(2 * nev + 1, 30))
        w = w[np.argsort(-np.abs(w))]
        return w[:nev]


def main():
    from models import model
    ap = argparse.ArgumentParser()
    ap.add_argument('model')
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--chis', type=int, nargs='+', default=[8, 12, 16, 24, 32, 48])
    ap.add_argument('--sweeps', type=int, default=6)
    ap.add_argument('--etol', type=float, default=1e-9)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--win', type=int, nargs='+', default=[32, 64, 128])
    ap.add_argument('--rmax', type=int, default=0)
    ap.add_argument('--out', default=None)
    ap.add_argument('--edgefield', action='store_true',
                    help='add f.sigma/2 on both end sites (full field on every site, as in d2_classification/dmrg.py)')
    a = ap.parse_args()
    h, r, c, note = model(a.model)
    W = mpo_from_term(h, r)
    edge = None
    if a.edgefield:
        import json as _j
        p = _j.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'candidates_params.json')))[a.model]
        edge = 0.5 * np.einsum('a,aij->ij', np.asarray(p['field'], float), PAULI)
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', f'fdmrg_{a.model}_L{a.L}_s{a.seed}.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    data = dict(model=a.model, note=note, c_exact=c, L=a.L, seed=a.seed, results=[])
    dm = FDMRG(W, a.L, seed=a.seed, edge=edge)
    data['edgefield'] = bool(a.edgefield)
    t0 = time.time()
    for chi in a.chis:
        Eprev = None
        hist = []
        for sw in range(a.sweeps * 4):
            E, tr = dm.sweep(chi)
            hist.append((E, tr))
            print(f"{a.model} L={a.L} chi={chi} sweep={sw} E={E:.12f} E/L={E/a.L:.10f} trunc={tr:.1e} "
                  f"Smid={dm.S[a.L//2-1]:.5f} t={time.time()-t0:.0f}s", flush=True)
            if sw + 1 >= a.sweeps and Eprev is not None and abs(E - Eprev) < a.etol * a.L:
                break
            Eprev = E
        rec = dict(chi=chi, E=E, sweeps=sw + 1, trunc=tr, S=dm.S.tolist(), hist=hist)
        mags = dm.local_ops(list(PAULI))
        rec['mag'] = mags.tolist()
        if r == 2:
            rec['Eb'] = dm.bond_energies(h).tolist()
        spec = {}
        for m in a.win:
            if m >= a.L - 4:
                continue
            j0 = a.L // 2 - m // 2
            w = dm.window_spectrum(j0, m)
            spec[str(m)] = [[float(z.real), float(z.imag)] for z in w]
            print(f"   window m={m}: |l1|={abs(w[1]):.6f} xi={-m/np.log(abs(w[1])/abs(w[0])):.3f} "
                  f"arg/pi={np.angle(w[1])/np.pi:.4f}", flush=True)
        rec['win'] = spec
        if a.rmax > 0:
            i0 = a.L // 2 - a.rmax // 2
            C = dm.correlations(i0, a.rmax, list(PAULI))
            rec['corr_i0'] = i0
            rec['corr'] = [[[C[x, y].real.tolist(), C[x, y].imag.tolist()] for y in range(3)] for x in range(3)]
        rec['time_s'] = time.time() - t0
        data['results'].append(rec)
        json.dump(data, open(out + '.tmp', 'w'))
        os.replace(out + '.tmp', out)


if __name__ == '__main__':
    main()

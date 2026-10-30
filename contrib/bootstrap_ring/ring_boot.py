"""Ring-aware bootstrap lower bounds on the ring ground energy E_0(L).

Parametrization. Real basis {P_mu} of d x d matrices, P_0 = 1, tr(P_mu^T P_nu) = d delta,
each symmetric or antisymmetric. The SDP variables are the expectation values
x_c = <P_mu1 (x) ... (x) P_mun> of operator strings, one variable per class of strings
under ring translation (sites mod L). Translation invariance, the r <-> L - r ring
identification and the compatibility of all marginals are therefore built in. Strings
with an odd number of antisymmetric factors have zero expectation in a real state; for a
real bond the twirled ground-space projector is real and translation invariant, so this
is a valid restriction.

Constraints (each is satisfied by the exact ring ground-space state; ed.py checks this):
  rho_k  >= 0       contiguous k-site marginal (k > L/2 makes it wrap around the ring)
  P_r    >= 0       two-block marginal on {0..s-1} u {r-1..r+s}, for every separation
                    whose sites do not fit in a k-arc (optional)
  G_loc  >= 0       <O_a^T [H,O_b]> for O on m_loc consecutive sites, with local
                    stationarity <[H,X]> = 0
  Gamma(q) >= 0     <O_q^dag O_q>/L, O_q = sum_j e^{iqj} O_j, O in the s-site algebra
  G(q)   >= 0       <O_q^dag [H,O_q]>/L, with its Hermiticity (bilocal stationarity)
at the momenta q = 2 pi m / L requested. For s >= 2 the set contains the energy density
modes h_q and their products with one-site operators.

Objective: minimise <h_01>. The optimum is a lower bound on E_0(L)/L.
"""
import os
for _n in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_n, '1')
import argparse, json, time
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
import cvxpy as cp

from models import bonds


def embed(op, pos, n, d):
    """Matrix on n sites acting as op on the ordered positions pos."""
    m = len(pos)
    rest = [i for i in range(n) if i not in pos]
    full = np.kron(op, np.eye(d ** (n - m)))
    perm = list(pos) + rest
    inv = list(np.argsort(perm))
    T = full.reshape((d,) * (2 * n)).transpose(inv + [n + i for i in inv])
    return T.reshape(d ** n, d ** n)


def site_basis(d):
    mats = [np.eye(d)]; par = [0]
    for j in range(1, d):
        v = np.zeros(d); v[:j] = 1; v[j] = -j
        mats.append(np.diag(v * np.sqrt(d / (j * (j + 1))))); par.append(0)
    for a in range(d):
        for b in range(a + 1, d):
            m = np.zeros((d, d)); m[a, b] = m[b, a] = np.sqrt(d / 2); mats.append(m); par.append(0)
            m = np.zeros((d, d)); m[a, b] = np.sqrt(d / 2); m[b, a] = -np.sqrt(d / 2); mats.append(m); par.append(1)
    return np.array(mats), np.array(par)


_SM = {}


def string_matrix(n, d):
    """Sparse (d^2n x d^2n): column mu = vec_C(P_mu1 (x) ... (x) P_mun)."""
    if (n, d) in _SM:
        return _SM[(n, d)]
    mats, _ = site_basis(d)
    Bm = sp.csr_matrix(mats.reshape(d * d, d * d).T)  # rows (a,b), cols mu
    M = Bm
    for _ in range(n - 1):
        M = sp.kron(M, Bm, format='csr')
    # rows are ordered (a1,b1,a2,b2,...); permute to (a1..an,b1..bn)
    D = d ** n
    order = [2 * i for i in range(n)] + [2 * i + 1 for i in range(n)]
    old = np.arange(D * D).reshape((d, d) * n).transpose(order).ravel()
    M = M[old, :]
    _SM[(n, d)] = M
    return M


class Ring:
    def __init__(self, h, d, L, k, s, m_loc, momenta, use_pair_psd=True, use_gram=True,
                 use_G=True, verbose=False, opset='full', pair_s=None):
        self.h, self.d, self.L, self.k, self.s = h, d, L, k, s
        self.ops = operator_set(h, d, s, opset)
        _, self.par = site_basis(d)
        self.reg = {}
        self.psd = []        # (name, list of (A_cols, A_rows, A_vals) deferred, ...)
        self.blocks = []     # deferred affine blocks: dict(kind, Cmat, U, shape)
        self.eq_rows = []    # deferred equality blocks (Cmat, U)
        self.log = verbose
        # Moment-matrix entries are string expectations, so they need no enclosing
        # density-matrix object; k only sets which strings are tied to rho_k >= 0.
        assert m_loc + 2 <= k
        # objects
        self.add_object('rho', list(range(k)))
        self.pairs = []
        ps = s if pair_s is None else pair_s
        self.pair_s = ps
        if s > 0 and use_pair_psd and (use_gram or use_G):
            for r in range(ps + 1, L - ps):
                sites = list(range(ps)) + [(r - 1 + j) % L for j in range(ps + 2)]
                if len(set(sites)) < len(sites) or self.arc_fits(sites):
                    continue
                self.add_object('P%d' % r, sites)
                self.pairs.append(r)
        # energy
        self.energy_block = self.functional(np.kron(np.eye(1), embed(self.h, [0, 1], 2, d))[None], [0, 1])
        self.local_ground(m_loc)
        self.momenta = momenta
        self.use_gram, self.use_G = use_gram, use_G
        self.g = {}
        if s > 0 and (use_gram or use_G):
            for kind, on in (('Gram', use_gram), ('G', use_G)):
                if on:
                    self.g[kind] = [self.entry_block(r, kind) for r in range(L)]
            if use_G:
                N = len(self.ops)
                for r in range(L // 2 + 1):
                    a = self.g['G'][r]; b = self.g['G'][(-r) % L]
                    # g_r[a,b] - g_{-r}[b,a] = 0
                    Tb = self.transpose_rows(b, N)
                    self.eq_rows.append(('diff', a, Tb))

    # ---------- strings and variables ----------
    def arc_fits(self, sites):
        k, L = self.k, self.L
        return any(all(((x + t) % L) < k for x in sites) for t in range(L))

    def cls(self, items):
        L = self.L
        best = None
        for p, _ in items:
            key = tuple(sorted((((x - p) % L), m) for x, m in items))
            if best is None or key < best:
                best = key
        return best

    def string_map(self, U, subset=None):
        """Map from string index (or position in subset) to variable index; odd strings
        and the identity map to nothing (the identity is handled as a constant)."""
        d, n = self.d, len(U)
        rows, cols = [], []
        idxs = np.arange(d ** (2 * n)) if subset is None else np.asarray(subset)
        if len(idxs) == 0:
            return (np.zeros(0, int), np.zeros(0, int), 0, None)
        mus = np.array(np.unravel_index(idxs, (d * d,) * n)).T
        ident = None
        for row, mu in enumerate(mus):
            if not mu.any():
                ident = row
                continue
            if self.par[mu].sum() % 2:
                continue
            items = [(U[i], int(mu[i])) for i in range(n) if mu[i]]
            c = self.cls(items)
            if c not in self.reg:
                self.reg[c] = len(self.reg)
            rows.append(row); cols.append(self.reg[c])
        return (np.array(rows, dtype=int), np.array(cols, dtype=int), len(idxs), ident)

    def functional(self, Ws, U):
        """Expectations of the operators Ws (batch, D, D) on ordered sites U, deferred.
        Only strings with a nonzero coefficient are kept."""
        d, n = self.d, len(U)
        SM = string_matrix(n, d)
        Wv = sp.csr_matrix(Ws.reshape(len(Ws), -1))
        C = sp.csr_matrix(Wv @ SM) / d ** n      # coefficient of each string
        C.data[np.abs(C.data) < 1e-13] = 0.0
        C.eliminate_zeros()
        used = np.unique(C.indices)
        C = sp.csr_matrix(C[:, used])
        return dict(C=C, S=self.string_map(U, used))

    def add_object(self, name, sites):
        d, n = self.d, len(sites)
        SM = string_matrix(n, d) / d ** n   # vec(rho) = SM @ <P_mu>
        self.blocks.append(dict(name=name, C=SM.tocsr(), S=self.string_map(sites), D=d ** n))

    def transpose_rows(self, blk, N):
        perm = np.arange(N * N).reshape(N, N).T.ravel()
        return dict(C=blk['C'][perm, :], S=blk['S'])

    def bond_sum(self, n, first, last):
        d = self.d
        H = np.zeros((d ** n, d ** n))
        for i in range(first, last):
            H += embed(self.h, [i, i + 1], n, d)
        return H

    def local_ground(self, m_loc):
        self.gloc = None
        if m_loc <= 0:
            return
        d = self.d
        n = m_loc + 2
        Hn = self.bond_sum(n, 0, n - 1)
        mid = list(range(1, m_loc + 1))
        D = d ** m_loc
        units = []
        for a in range(D):
            for b in range(D):
                e = np.zeros((D, D)); e[a, b] = 1; units.append(embed(e, mid, n, d))
        comm = [Hn @ o - o @ Hn for o in units]
        Ws = np.array([oa.T @ cb for oa in units for cb in comm])
        self.gloc = self.functional(Ws, list(range(n)))
        self.gloc_N = len(units)
        # stationarity <[H, X]> = 0 for every X on the middle sites
        self.eq_rows.append(('zero', self.functional(np.array(comm), list(range(n)))))

    def entry_block(self, r, kind):
        d, s, L = self.d, self.s, self.L
        A = list(range(s))
        B = [(r - 1 + j) % L for j in range(s + 2)] if kind == 'G' else [(r + j) % L for j in range(s)]
        U = []
        for x in A + B:
            if x not in U:
                U.append(x)
        n = len(U)
        posA = [U.index(x) for x in A]; posB = [U.index(x) for x in B]
        units = self.ops
        if kind == 'G':
            HB = self.bond_sum(s + 2, 0, s + 1)
            mid = list(range(1, s + 1))
            Bops = []
            for e in units:
                o = embed(e, mid, s + 2, d)
                Bops.append(embed(HB @ o - o @ HB, posB, n, d))
        else:
            Bops = [embed(e, posB, n, d) for e in units]
        Aops = [embed(e.T, posA, n, d) for e in units]
        Ws = np.array([oa @ ob for oa in Aops for ob in Bops])
        return self.functional(Ws, U)

    # ---------- assembly ----------
    def affine(self, blk):
        """Return sparse A (rows x nvar) and constant b for a deferred block."""
        rows, cols, nS, ident = blk['S']
        S = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(nS, self.nvar))
        A = blk['C'] @ S
        if ident is None:
            b = np.zeros(blk['C'].shape[0])
        else:
            b = np.asarray(blk['C'][:, ident].todense()).ravel()
        return sp.csr_matrix(A), b

    def build(self):
        self.nvar = len(self.reg)
        x = cp.Variable(self.nvar)
        self.x = x
        cons = []
        used = np.zeros(self.nvar, bool)
        psd = []                      # (A, b, D) symmetric affine maps, vec_C order

        def clean(A, b):
            A = sp.csr_matrix(A); A.sum_duplicates()
            A.data[np.abs(A.data) < 1e-12] = 0.0; A.eliminate_zeros()
            b = np.where(np.abs(b) < 1e-12, 0.0, b)
            return A, b

        def sym(A, b, D, sign=1):
            T = np.arange(D * D).reshape(D, D).T.ravel()
            return clean((A + sign * A[T]) / 2, (b + sign * b[T]) / 2)

        for blk in self.blocks:
            A, b = self.affine(blk)
            psd.append(sym(A, b, blk['D']) + (blk['D'],))
        A, b = self.affine(self.energy_block)
        self.energy = A[0] @ x + b[0]
        used[A[0].indices] = True
        if self.gloc is not None:
            A, b = self.affine(self.gloc)
            N = self.gloc_N
            self.Gloc_expr = cp.reshape(A @ x + b, (N, N), order='C')
            psd.append(sym(A, b, N) + (N,))
            # antisymmetric part of G_loc is local stationarity of O_a^T O_b
            Aa, ba = sym(A, b, N, -1)
            self.eq_rows.append(('raw', Aa, ba))
        # momentum-space moment matrices, as real 2N x 2N embeddings of Hermitian parts
        self.mom = {}
        for kind, blks in self.g.items():
            N = len(self.ops)
            aff = [self.affine(bk) for bk in blks]
            for m in self.momenta:
                q = 2 * np.pi * m / self.L
                cr = [np.cos(q * r) for r in range(self.L)]
                sr = [np.sin(q * r) for r in range(self.L)]
                cr = [0.0 if abs(c) < 1e-12 else c for c in cr]
                sr = [0.0 if abs(c) < 1e-12 else c for c in sr]
                ARe = sum(c * a[0] for c, a in zip(cr, aff)); bRe = sum(c * a[1] for c, a in zip(cr, aff))
                AIm = sum(c * a[0] for c, a in zip(sr, aff)); bIm = sum(c * a[1] for c, a in zip(sr, aff))
                self.mom[(kind, m)] = (cp.reshape(ARe @ x + bRe, (N, N), order='C'),
                                       cp.reshape(AIm @ x + bIm, (N, N), order='C'))
                Rs, bR = sym(ARe, bRe, N)
                Ia, bI = sym(AIm, bIm, N, -1)
                # E = [[Rs, -Ia], [Ia, Rs]] in vec_C order of a 2N x 2N matrix
                idx = np.arange(4 * N * N).reshape(2 * N, 2 * N)
                rows_R1 = idx[:N, :N].ravel(); rows_R2 = idx[N:, N:].ravel()
                rows_I1 = idx[N:, :N].ravel(); rows_I2 = idx[:N, N:].ravel()
                P = lambda rows: sp.csr_matrix((np.ones(N * N), (rows, np.arange(N * N))), shape=(4 * N * N, N * N))
                AE = P(rows_R1) @ Rs + P(rows_R2) @ Rs + P(rows_I1) @ Ia - P(rows_I2) @ Ia
                bE = P(rows_R1) @ bR + P(rows_R2) @ bR + P(rows_I1) @ bI - P(rows_I2) @ bI
                psd.append(clean(AE, bE) + (2 * N,))
        for A, b, D in psd:
            used[A.indices] = True
            M = cp.reshape(A @ x + b, (D, D), order='C')
            cons.append((M + M.T) / 2 >> 0)
        # equalities
        Arows, brows = [], []
        for item in self.eq_rows:
            if item[0] == 'zero':
                A, b = self.affine(item[1])
            elif item[0] == 'diff':
                A1, b1 = self.affine(item[1]); A2, b2 = self.affine(item[2])
                A, b = A1 - A2, b1 - b2
            else:
                A, b = item[1], item[2]
            Arows.append(sp.csr_matrix(A)); brows.append(b)
        self.n_eq_raw = 0
        if Arows:
            A = sp.vstack(Arows).tocsr(); b = np.concatenate(brows)
            A, b = independent_rows(A, b)
            self.n_eq_raw = sum(a.shape[0] for a in Arows)
            self.n_eq = A.shape[0]
            # string variables that enter no PSD block and no objective are fixed by the
            # equalities alone or are irrelevant; pin the irrelevant ones to zero so the
            # KKT system is nonsingular (this does not change the optimum).
            if A.shape[0]:
                cons.append(A @ x + b == 0)
                used[A.indices] = True
        free = np.where(~used)[0]
        self.n_unused = len(free)
        if len(free):
            cons.append(x[free] == 0)
        self.cons = cons
        return cons

    def solve(self, solver='CLARABEL'):
        t0 = time.time()
        self.build()
        tb = time.time() - t0
        t = time.time()
        attempts = [(solver, 0.0)]
        if solver == 'CLARABEL':
            # A tiny ridge term makes the KKT system nonsingular when some string
            # variables enter the constraints only in fixed combinations. The reported
            # value is the energy at the returned point; the ridge shifts the optimum by
            # at most reg * |x_opt|^2.
            attempts = [('CLARABEL', 0.0), ('CLARABEL', 1e-9), ('SCS', 0.0)]
        else:
            attempts += [('SCS', 0.0)]
        for used, reg in attempts:
            obj = self.energy if reg == 0 else self.energy + reg * cp.sum_squares(self.x)
            prob = cp.Problem(cp.Minimize(obj), self.cons)
            try:
                if used == 'SCS':
                    prob.solve(solver='SCS', eps_abs=1e-9, eps_rel=1e-9, max_iters=100000, verbose=self.log)
                else:
                    # stronger static regularization is needed when the string
                    # parametrization leaves near-degenerate directions (s = 2 runs)
                    prob.solve(solver=used, verbose=self.log, static_regularization_constant=1e-7)
                if prob.status in ('optimal', 'optimal_inaccurate'):
                    break
            except cp.error.SolverError:
                continue
        self.reg_used = reg
        val = float(np.ravel(self.energy.value)[0]) if self.x.value is not None else None
        return prob, dict(status=prob.status, value=val, solver=used, ridge=reg,
                          x_norm2=float(np.sum(self.x.value ** 2)) if self.x.value is not None else None,
                          seconds=time.time() - t, assemble_seconds=tb, nvar=self.nvar,
                          n_eq=getattr(self, 'n_eq', 0), n_eq_raw=self.n_eq_raw)


def operator_set(h, d, s, opset):
    """Local operators O (on s sites) whose momentum modes O_q enter Gamma(q) and G(q)."""
    D = d ** s
    units = []
    for a in range(D):
        for b in range(D):
            e = np.zeros((D, D)); e[a, b] = 1; units.append(e)
    if opset == 'full':
        return units
    if opset == 'onesite_h':
        assert s == 2
        one = []
        for a in range(d):
            for b in range(d):
                e = np.zeros((d, d)); e[a, b] = 1; one.append(np.kron(e, np.eye(d)))
        return one + [h.copy()]
    if opset == 'h':
        assert s == 2
        return [np.eye(D), h.copy()]
    raise ValueError(opset)


def independent_rows(A, b, tol=1e-7):
    """Remove zero and linearly dependent rows of A x + b = 0 (dense QR on used columns)."""
    A = sp.csr_matrix(A)
    A.sum_duplicates()
    A.data[np.abs(A.data) < 1e-12] = 0.0          # floating cancellation residue
    b = np.where(np.abs(b) < 1e-12, 0.0, b)
    A.eliminate_zeros(); A.sort_indices()
    nz = np.diff(A.indptr) > 0
    assert np.all(b[~nz] == 0), 'inconsistent constant row'
    A = A[nz]; b = b[nz]
    if A.shape[0] == 0:
        return A, b
    used = np.unique(A.indices)
    Ad = np.asarray(A[:, used].todense())
    Ab = np.hstack([Ad, b[:, None]])
    scale = np.linalg.norm(Ab, axis=1)
    Ab = Ab / scale[:, None]
    Q, R, piv = sla.qr(Ab.T, mode='economic', pivoting=True)
    diag = np.abs(np.diag(R))
    rank = int(np.sum(diag > tol * diag[0]))
    keep = np.sort(piv[:rank])
    Ar = sp.csr_matrix(sp.diags(1 / scale[keep]) @ A[keep])
    br = b[keep] / scale[keep]
    return Ar, br


def momenta_list(L, mode):
    if mode == 'all':
        return list(range(L // 2 + 1))
    if mode == 'none':
        return []
    return [int(x) for x in mode.split(',')]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--model', required=True)
    p.add_argument('--L', type=int, required=True)
    p.add_argument('--k', type=int, required=True)
    p.add_argument('--s', type=int, default=1)
    p.add_argument('--mloc', type=int, default=None)
    p.add_argument('--momenta', default='all')
    p.add_argument('--no-pair-psd', action='store_true')
    p.add_argument('--no-gram', action='store_true')
    p.add_argument('--no-G', action='store_true')
    p.add_argument('--solver', default='CLARABEL')
    p.add_argument('--opset', default='full', choices=['full', 'onesite_h', 'h'])
    p.add_argument('--pair-s', type=int, default=None, help='block size of the two-block PSD objects (default s)')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    if 'SLURM_JOB_ID' not in os.environ:
        p.error('run inside a SLURM allocation')
    B = {nm: (d, h) for nm, d, h in bonds()}
    d, h = B[a.model]
    mloc = a.mloc if a.mloc is not None else a.k - 2
    s = a.s if not (a.no_gram and a.no_G) else 0
    t = time.time()
    R = Ring(h, d, a.L, a.k, s, mloc, momenta_list(a.L, a.momenta),
             use_pair_psd=not a.no_pair_psd, use_gram=not a.no_gram, use_G=not a.no_G,
             opset=a.opset, pair_s=a.pair_s)
    build = time.time() - t
    prob, info = R.solve(a.solver)
    row = dict(model=a.model, d=d, L=a.L, k=a.k, s=a.s, mloc=mloc, momenta=a.momenta, opset=a.opset, pair_s=R.pair_s,
               pair_psd=not a.no_pair_psd, gram=not a.no_gram, G=not a.no_G,
               n_pairs=len(R.pairs), build_seconds=build, **info)
    if info['value'] is not None:
        row['bound_per_site'] = info['value']
        row['bound_total'] = a.L * info['value']
    with open(a.out, 'w') as f:
        json.dump(row, f)
    print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()

"""C137 window lemma: shared code.

Conventions follow C135 core.py: S = sigma/2, XXZ bond = SxSx + SySy + Delta SzSz,
TFIM bond = -XX - g (ZI + IZ)/2 in Pauli matrices. Basis ordering: site 1 is the
most significant bit.

A window on n sites is W(c, u) = sum_{j=1}^{n-1} c_j h_{j,j+1} + sum_{k=1}^n u_k,
with sum_j c_j = 1 and sum_k u_k = 0 (u_k traceless one-site operators).
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
import json, time
from pathlib import Path
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla

HERE = Path(__file__).resolve().parent
C135 = HERE.parent / 'work_C135'

I2 = np.eye(2)
PX = np.array([[0, 1], [1, 0]], complex)
PY = np.array([[0, -1j], [1j, 0]], complex)
PZ = np.array([[1, 0], [0, -1]], complex)
PAULI = [PX, PY, PZ]


def kron(*a):
    out = np.array([[1.0 + 0j]])
    for x in a:
        out = np.kron(out, x)
    return out


def xxz(delta, b=0.0, dm=0.0):
    S = [p / 2 for p in PAULI]
    h = kron(S[0], S[0]) + kron(S[1], S[1]) + delta * kron(S[2], S[2])
    h = h + dm * (kron(S[0], S[1]) - kron(S[1], S[0]))
    h = h + b * (kron(S[2], I2) + kron(I2, S[2])) / 2
    return h


def xyz(jx, jy, jz):
    S = [p / 2 for p in PAULI]
    return jx * kron(S[0], S[0]) + jy * kron(S[1], S[1]) + jz * kron(S[2], S[2])


def tfim(g, lam=0.0):
    """-XX - g (ZI+IZ)/2 + lam ZZ (lam = 0 is the free TFIM)."""
    return -kron(PX, PX) - g * (kron(PZ, I2) + kron(I2, PZ)) / 2 + lam * kron(PZ, PZ)


def herm(a):
    return (a + a.conj().T) / 2


def gauge_M(h):
    """Gauge-optimized half-spread M(h) by a small SDP, with a certified lower bound.
    Returns dict(M_up, M_low, o) where o are the one-site gauge Pauli coefficients."""
    import cvxpy as cp
    G = [kron(p, I2) - kron(I2, p) for p in PAULI]
    o = cp.Variable(3); a = cp.Variable(); m = cp.Variable()
    q = h + a * np.eye(4) + sum(o[j] * G[j] for j in range(3))
    cons = [m * np.eye(4) - q >> 0, m * np.eye(4) + q >> 0]
    p = cp.Problem(cp.Minimize(m), cons)
    try:
        p.solve(solver='CLARABEL', tol_gap_abs=1e-11, tol_feas=1e-11, tol_gap_rel=1e-11)
    except Exception:
        p.solve(solver='SCS', eps=1e-10, max_iters=50000)
    rep = h + sum(o.value[j] * G[j] for j in range(3))
    ev = la.eigvalsh(herm(rep)); up = (ev[-1] - ev[0]) / 2
    D = herm(cons[0].dual_value - cons[1].dual_value)
    D -= np.trace(D) * np.eye(4) / 4
    for g in G:
        D -= np.trace(D @ g) / np.trace(g @ g) * g
    D = herm(D)
    D /= max(1.0, float(np.abs(la.eigvalsh(D)).sum()))
    low = float(np.trace(h @ D).real)
    return dict(M_up=float(up), M_low=low, o=[float(x) for x in o.value])


# ---------------------------------------------------------------- sparse window operators
def site_op(op, k, n):
    """Sparse operator op (2x2) on site k (0-based) of n sites."""
    return sp.kron(sp.kron(sp.identity(2 ** k, format='csr'), sp.csr_matrix(op)),
                   sp.identity(2 ** (n - k - 1), format='csr'), format='csr')


def bond_op(h, j, n):
    """Sparse operator h (4x4) on sites j, j+1 (0-based)."""
    return sp.kron(sp.kron(sp.identity(2 ** j, format='csr'), sp.csr_matrix(h)),
                   sp.identity(2 ** (n - j - 2), format='csr'), format='csr')


class WindowFamily:
    """Affine family A(z) = A0 + sum_i z_i A_i of window operators on n sites.

    Parametrization: c = cbar + Bc z_c (sum zero directions), u_k = sum_a y_{k,a} sigma^a with
    y = By z_y (sum over k zero for each a). Optional reflection symmetry (c_j = c_{n-j},
    y_k = y_{n+1-k}) and a restricted set of Pauli directions for u.
    """

    def __init__(self, h, n, udirs=(0, 1, 2), reflect=False, cfix=None):
        self.h = np.asarray(h, complex); self.n = n
        self.real = np.max(np.abs(self.h.imag)) < 1e-14 and all(a != 1 for a in udirs)
        dt = float if self.real else complex
        hh = self.h.real if self.real else self.h
        self.bonds = [bond_op(hh, j, n).astype(dt) for j in range(n - 1)]
        self.udirs = list(udirs)
        self.sites = {a: [site_op(PAULI[a].real if self.real else PAULI[a], k, n).astype(dt)
                          for k in range(n)] for a in self.udirs}
        m = n - 1
        # c directions
        if cfix is not None:
            self.cbar = np.asarray(cfix, float); Bc = np.zeros((m, 0))
        else:
            self.cbar = np.ones(m) / m
            Bc = sum_zero_basis(m, reflect)
        # u directions, per Pauli direction a
        By = sum_zero_basis(n, reflect)
        self.Bc = Bc; self.By = By
        self.nz = Bc.shape[1] + len(self.udirs) * By.shape[1]
        self.dt = dt

    def split(self, z):
        nc = self.Bc.shape[1]
        c = self.cbar + self.Bc @ z[:nc]
        ys = {}
        rest = z[nc:]
        nb = self.By.shape[1]
        for i, a in enumerate(self.udirs):
            ys[a] = self.By @ rest[i * nb:(i + 1) * nb]
        return c, ys

    def matrix(self, z):
        c, ys = self.split(z)
        A = sum(cj * B for cj, B in zip(c, self.bonds))
        for a, y in ys.items():
            for k in range(self.n):
                if abs(y[k]) > 0:
                    A = A + y[k] * self.sites[a][k]
        return A.tocsr()

    def direction_ops(self):
        """List of sparse A_i with A(z) = A(0) + sum z_i A_i."""
        ops = []
        for col in self.Bc.T:
            ops.append(sum(cj * B for cj, B in zip(col, self.bonds) if abs(cj) > 0).tocsr())
        for a in self.udirs:
            for col in self.By.T:
                ops.append(sum(yk * S for yk, S in zip(col, self.sites[a]) if abs(yk) > 0).tocsr())
        return ops


def sum_zero_basis(m, reflect):
    """Orthonormal basis of vectors in R^m with zero sum (and reflection symmetry if asked)."""
    if m <= 1:
        return np.zeros((m, 0))
    if reflect:
        cols = []
        for j in range((m + 1) // 2):
            v = np.zeros(m); v[j] += 1; v[m - 1 - j] += 1
            cols.append(v)
        V = np.array(cols).T
    else:
        V = np.eye(m)
    V = V - V.mean(axis=0, keepdims=True)
    U, s, _ = la.svd(V, full_matrices=False)
    return U[:, s > 1e-10]


def lowest(A, k=4, v0=None, tol=1e-12):
    N = A.shape[0]
    if N <= 1024:
        ev, V = la.eigh(A.toarray())
        return ev[:k], V[:, :k]
    k = min(k, N - 2)
    ev, V = sla.eigsh(A, k=k, which='SA', tol=tol, v0=v0, ncv=max(2 * k + 8, 24))
    ix = np.argsort(ev)
    return ev[ix], V[:, ix]


def optimize_window(fam, z0=None, maxit=200, tol=1e-10, k=4, rmax=48, verbose=False, mu0=None):
    """Maximize lambda_min(A(z)) (concave) by a proximal spectral-bundle method.

    Returns dict with best value, z, c, y, number of iterations and an a-posteriori
    model gap (max of the subspace model in a box around the best point minus the value)."""
    import cvxpy as cp
    p = fam.nz
    z = np.zeros(p) if z0 is None else np.array(z0, float)
    A = fam.matrix(z)
    ev, V = lowest(A, k)
    f = ev[0]
    if p == 0:
        c, ys = fam.split(z)
        return dict(value=float(f), z=z, c=c, y=ys, iters=0, pred_gap=0.0, evals=ev.tolist(), hist=[])
    ops = fam.direction_ops()
    Q = V.copy()
    mu = mu0 if mu0 is not None else 1.0
    hist = []
    it = 0
    pred = np.inf
    v0 = V[:, 0]
    for it in range(maxit):
        Q, _ = la.qr(Q, mode='economic')
        r = Q.shape[1]
        A0 = Q.conj().T @ (A @ Q)
        Ai = [Q.conj().T @ (Op @ Q) for Op in ops]
        # model in coordinates relative to current z: A(z + d) = A(z) + sum d_i A_i
        d = cp.Variable(p); t = cp.Variable()
        M = (A0 + sum(d[i] * Ai[i] for i in range(p))) if False else None
        # real embedding of Hermitian LMI
        def emb(X):
            X = herm(X)
            if np.iscomplexobj(X) and np.max(np.abs(X.imag)) > 1e-15:
                return np.block([[X.real, -X.imag], [X.imag, X.real]])
            return np.block([[X.real, np.zeros_like(X.real)], [np.zeros_like(X.real), X.real]])
        E0 = emb(A0); Ei = [emb(X) for X in Ai]
        expr = E0 + sum(d[i] * Ei[i] for i in range(p)) - t * np.eye(2 * r)
        prob = cp.Problem(cp.Maximize(t - mu / 2 * cp.sum_squares(d)), [(expr + expr.T) / 2 >> 0])
        ok = False
        for solver, opts in [('CLARABEL', {}), ('SCS', dict(eps=1e-9, max_iters=50000)), ('CVXOPT', {})]:
            try:
                prob.solve(solver=solver, **opts)
                if t.value is not None and d.value is not None:
                    ok = True; break
            except Exception:
                pass
        if not ok:
            mu = min(mu * 4, 1e8)
            if mu >= 1e8:
                break
            continue
        dval = np.array(d.value).ravel(); mval = float(t.value)
        pred = mval - f
        if pred < tol * max(1, abs(f)):
            break
        znew = z + dval
        Anew = fam.matrix(znew)
        evn, Vn = lowest(Anew, k, v0=v0)
        fn = evn[0]
        if fn - f >= 0.2 * pred:
            z, A, f, ev, V = znew, Anew, fn, evn, Vn
            v0 = V[:, 0]
            mu = max(mu / 2, 1e-6)
        else:
            mu = min(mu * 4, 1e8)
        Q = np.hstack([V, Vn, Q])[:, :rmax]
        hist.append((it, float(f), float(pred), mu, r))
        if verbose:
            print(it, f, pred, mu, r, flush=True)
    c, ys = fam.split(z)
    return dict(value=float(f), z=z, c=c, y={a: ys[a] for a in ys}, iters=it, pred_gap=float(pred),
                evals=ev.tolist(), hist=hist)


def window_value(h, c, n, udirs=(0, 1, 2), reflect=False, optimize_u=True, k=4):
    """lambda_min for fixed c, with u optimized (concave in u) or u = 0."""
    fam = WindowFamily(h, n, udirs=udirs if optimize_u else (), reflect=reflect, cfix=c)
    if not optimize_u or fam.nz == 0:
        ev, _ = lowest(fam.matrix(np.zeros(fam.nz)), k)
        return dict(value=float(ev[0]))
    return optimize_window(fam, k=k)


# ---------------------------------------------------------------- standard windows
def ssd_window(n):
    j = np.arange(1, n)
    c = np.sin(np.pi * j / n) ** 2
    return c / c.sum()


def ssd_site_weights(n):
    k = np.arange(1, n + 1)
    w = np.sin(np.pi * (k - 0.5) / n) ** 2
    return w / w.sum()


def sine_window(n):
    """C135 window a_j = v_j v_{j+1}/s, v_j = sqrt(2/(n+1)) sin(pi j/(n+1))."""
    v = np.sqrt(2 / (n + 1)) * np.sin(np.pi * np.arange(1, n + 1) / (n + 1))
    a = v[:-1] * v[1:]
    return a / a.sum()


def sine_site_weights(n):
    v = np.sqrt(2 / (n + 1)) * np.sin(np.pi * np.arange(1, n + 1) / (n + 1))
    return v ** 2


def write_json(path, x):
    Path(path).write_text(json.dumps(x, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n')

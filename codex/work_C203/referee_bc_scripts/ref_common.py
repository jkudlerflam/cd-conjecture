"""Referee C203: independent ED helpers (not reusing C201 scripts)."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def digits(L, d):
    idx = np.arange(d ** L)
    D = np.empty((d ** L, L), dtype=np.int64)
    x = idx.copy()
    for k in range(L - 1, -1, -1):
        D[:, k] = x % d
        x //= d
    return D  # D[state, site], site 0 is most significant


def two_site_op(h, i, j, L, d, D=None):
    """Sparse operator of 2-site h (d^2 x d^2, first factor = site i) acting on sites i, j."""
    if D is None:
        D = digits(L, d)
    N = d ** L
    pw = d ** (L - 1 - np.arange(L))
    rows, cols, vals = [], [], []
    a = D[:, i]
    b = D[:, j]
    col_ab = a * d + b
    base = np.arange(N) - a * pw[i] - b * pw[j]
    for out in range(d * d):
        ap, bp = divmod(out, d)
        coeff = h[out, col_ab]
        nz = np.abs(coeff) > 0
        if not np.any(nz):
            continue
        rows.append(base[nz] + ap * pw[i] + bp * pw[j])
        cols.append(np.arange(N)[nz])
        vals.append(coeff[nz])
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))


def ring_H(h, L, d):
    D = digits(L, d)
    H = sp.csr_matrix((d ** L, d ** L), dtype=complex)
    for j in range(L):
        H = H + two_site_op(h, j, (j + 1) % L, L, d, D)
    return H


def lowest(H, k=1):
    n = H.shape[0]
    if n <= 3000:
        w = np.linalg.eigvalsh(H.toarray() if sp.issparse(H) else H)
        return w[:k]
    w = spla.eigsh(H, k=k, which='SA', tol=1e-12, maxiter=100000)[0]
    return np.sort(w)


def gauge_M(h, d):
    """min over Hermitian X of half spectral spread of h + X(x)1 - 1(x)X, via SDP."""
    import cvxpy as cp
    X = cp.Variable((d, d), hermitian=True)
    s = cp.Variable()
    r = cp.Variable()
    I = np.eye(d)
    G = h + cp.kron(X, I) - cp.kron(I, X)
    G = (G + G.H) / 2
    cons = [G << s * np.eye(d * d), G >> r * np.eye(d * d)]
    prob = cp.Problem(cp.Minimize(s - r), cons)
    prob.solve(solver=cp.CLARABEL)
    return (s.value - r.value) / 2


def haar_unitary(n, rng):
    Z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))

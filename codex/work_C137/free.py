"""Exact window energies for quadratic qubit bonds by the Majorana method, to large n.

Majoranas a_j = (prod_{l<j} Z_l) X_j, b_j = (prod_{l<j} Z_l) Y_j on an open chain. Then
  Z_j = -i a_j b_j,  X_j X_{j+1} = -i b_j a_{j+1},  Y_j Y_{j+1} = i a_j b_{j+1},
  X_j Y_{j+1} = -i b_j b_{j+1},  Y_j X_{j+1} = i a_j a_{j+1}.
A quadratic H = const + (i/4) gamma^T A gamma, A real antisymmetric, has ground energy
const - (1/4) ||A||_* (nuclear norm). Maximizing lambda_min over windows is therefore a
nuclear-norm minimization, a convex problem.
"""
import numpy as np
import scipy.linalg as la
from scipy.integrate import quad
from common import PAULI, I2, kron

LAB = {'I': I2, 'X': PAULI[0], 'Y': PAULI[1], 'Z': PAULI[2]}


def pauli_coeffs(h):
    out = {}
    for a in 'IXYZ':
        for b in 'IXYZ':
            out[a + b] = np.trace(kron(LAB[a], LAB[b]) @ h).real / 4
    return out


def is_quadratic(h, tol=1e-12):
    p = pauli_coeffs(h)
    bad = [k for k in p if k not in ('II', 'ZI', 'IZ', 'XX', 'XY', 'YX', 'YY') and abs(p[k]) > tol]
    return len(bad) == 0, bad


# antisymmetric two-index contributions: term coefficient t times operator = t * (i * s) gamma_p gamma_q
# with (i/4)(A_pq gamma_p gamma_q + A_qp gamma_q gamma_p) = (i/2) A_pq gamma_p gamma_q, so A_pq += 2 t s.
def build_A(n, bond_coef, site_z):
    """bond_coef: list of dicts {XX,XY,YX,YY} per bond j (0..n-2); site_z: array of Z coefficients."""
    A = np.zeros((2 * n, 2 * n))
    a = lambda j: 2 * j
    b = lambda j: 2 * j + 1
    def add(p, q, val):
        A[p, q] += val; A[q, p] -= val
    for j, bc in enumerate(bond_coef):
        add(b(j), a(j + 1), 2 * bc['XX'] * (-1))
        add(a(j), b(j + 1), 2 * bc['YY'] * (+1))
        add(b(j), b(j + 1), 2 * bc['XY'] * (-1))
        add(a(j), a(j + 1), 2 * bc['YX'] * (+1))
    for k, z in enumerate(site_z):
        add(a(k), b(k), 2 * z * (-1))
    return A


def window_energy(h, c, site_extra=None, onsite_weights=None):
    """lambda_min of sum_j c_j h_{j,j+1} + sum_k site_extra_k Z_k.
    If onsite_weights is given, the Z part of the bond is instead placed with total weight
    onsite_weights_k * (z1 + z2) on site k (this differs from the bond-induced weights by a
    sum-zero single-site term when sum(onsite_weights) = sum(c))."""
    p = pauli_coeffs(h)
    n = len(c) + 1
    bonds = [{k: cj * p[k] for k in ('XX', 'XY', 'YX', 'YY')} for cj in c]
    z = np.zeros(n)
    if onsite_weights is None:
        for j, cj in enumerate(c):
            z[j] += cj * p['ZI']; z[j + 1] += cj * p['IZ']
    else:
        z += np.asarray(onsite_weights) * (p['ZI'] + p['IZ'])
    if site_extra is not None:
        z += site_extra
    A = build_A(n, bonds, z)
    const = p['II'] * np.sum(c)
    return const - 0.25 * np.sum(np.abs(la.eigvalsh(1j * A)))


def e0_quadratic(h):
    """Energy density on Z from the 2x2 Majorana symbol."""
    p = pauli_coeffs(h)
    # infinite chain: A_0 (on-site 2x2 block for (a_j,b_j)), A_1 block coupling site j to j+1
    zt = p['ZI'] + p['IZ']
    A0 = np.array([[0, -2 * zt], [2 * zt, 0]])
    A1 = np.zeros((2, 2))  # rows (a_j, b_j), cols (a_{j+1}, b_{j+1})
    A1[1, 0] += -2 * p['XX']; A1[0, 1] += 2 * p['YY']; A1[1, 1] += -2 * p['XY']; A1[0, 0] += 2 * p['YX']
    def f(k):
        # (i/4) sum_{j,l} gamma_j^T A_{j-l} gamma_l ; symbol Ahat(k) = A0 + A1 e^{ik} - A1^T e^{-ik}
        Ah = A0 + A1 * np.exp(1j * k) - A1.T * np.exp(-1j * k)
        return np.sum(np.abs(la.eigvalsh(1j * Ah)))
    val = quad(f, -np.pi, np.pi, limit=400, epsabs=1e-14, epsrel=1e-13)[0] / (2 * np.pi)
    return p['II'] - 0.25 * val


def optimize_free(h, n, cfix=None, reflect=True, solver='CLARABEL'):
    """Maximize lambda_min over windows (c, y Z-fields) for a quadratic bond whose one-site gauge
    is along Z (true for the TFIM and U(1) bonds). Nuclear-norm minimization by cvxpy."""
    import cvxpy as cp
    p = pauli_coeffs(h)
    m = n - 1
    if cfix is None:
        c = cp.Variable(m)
        cons = [cp.sum(c) == 1]
        if reflect:
            cons += [c[j] == c[m - 1 - j] for j in range(m // 2)]
    else:
        c = np.asarray(cfix, float); cons = []
    y = cp.Variable(n)
    cons += [cp.sum(y) == 0]
    if reflect:
        cons += [y[k] == y[n - 1 - k] for k in range(n // 2)]
    # assemble A linearly
    E = lambda p_, q_: _E(2 * n, p_, q_)
    terms = []
    a = lambda j: 2 * j
    b = lambda j: 2 * j + 1
    Aexpr = 0
    for j in range(m):
        cj = c[j]
        Bm = (-2 * p['XX']) * E(b(j), a(j + 1)) + (2 * p['YY']) * E(a(j), b(j + 1)) \
            + (-2 * p['XY']) * E(b(j), b(j + 1)) + (2 * p['YX']) * E(a(j), a(j + 1)) \
            + (-2 * p['ZI']) * E(a(j), b(j)) + (-2 * p['IZ']) * E(a(j + 1), b(j + 1))
        Aexpr = Aexpr + cj * Bm
    for k in range(n):
        Aexpr = Aexpr + y[k] * (-2) * E(a(k), b(k))
    # nuclear norm of i A equals that of A (A real antisymmetric: singular values = |eigs of iA|)
    prob = cp.Problem(cp.Minimize(cp.normNuc(Aexpr)), cons)
    prob.solve(solver=solver)
    cval = c.value if cfix is None else c
    lam = p['II'] - 0.25 * prob.value
    # recompute exactly
    lam2 = window_energy(h, cval, site_extra=y.value)
    return dict(value=float(lam2), solver_value=float(lam), c=np.array(cval), y=np.array(y.value), status=prob.status)


def _E(N, p, q):
    M = np.zeros((N, N)); M[p, q] = 1; M[q, p] = -1
    return M


def optimize_free_lbfgs(h, n, cfix=None, reflect=True, eps_list=(1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 0.0)):
    """Same problem as optimize_free, by L-BFGS on the smoothed nuclear norm sum_i sqrt(mu_i^2 + eps^2)
    (mu_i the eigenvalues of i A), with eps decreased to 0. Fast for n in the hundreds."""
    from scipy.optimize import minimize
    from common import sum_zero_basis
    p = pauli_coeffs(h)
    m = n - 1
    a = lambda j: 2 * j
    b = lambda j: 2 * j + 1
    # elementary matrices: bond j operator Bm_j, site k operator Sk
    def bond_mat(j):
        M = np.zeros((2 * n, 2 * n))
        for (pp, qq, val) in [(b(j), a(j + 1), -2 * p['XX']), (a(j), b(j + 1), 2 * p['YY']),
                              (b(j), b(j + 1), -2 * p['XY']), (a(j), a(j + 1), 2 * p['YX']),
                              (a(j), b(j), -2 * p['ZI']), (a(j + 1), b(j + 1), -2 * p['IZ'])]:
            M[pp, qq] += val; M[qq, pp] -= val
        return M
    def site_mat(k):
        M = np.zeros((2 * n, 2 * n)); M[a(k), b(k)] = -2; M[b(k), a(k)] = 2
        return M
    Bms = [bond_mat(j) for j in range(m)]
    Sks = [site_mat(k) for k in range(n)]
    if cfix is None:
        cbar = np.ones(m) / m; Bc = sum_zero_basis(m, reflect)
    else:
        cbar = np.asarray(cfix, float); Bc = np.zeros((m, 0))
    By = sum_zero_basis(n, reflect)
    A0 = sum(cj * Bm for cj, Bm in zip(cbar, Bms))
    dirs = [sum(col[j] * Bms[j] for j in range(m)) for col in Bc.T] + [sum(col[k] * Sks[k] for k in range(n)) for col in By.T]
    P = len(dirs)
    D = np.array([d_.ravel() for d_ in dirs]) if P else np.zeros((0, 4 * n * n))
    def Aof(z):
        return A0 + (z @ D).reshape(2 * n, 2 * n) if P else A0
    def f(z, eps):
        A = Aof(z)
        mu, V = la.eigh(1j * A)
        if eps > 0:
            w = mu / np.sqrt(mu ** 2 + eps ** 2); val = np.sum(np.sqrt(mu ** 2 + eps ** 2))
        else:
            w = np.sign(mu); val = np.sum(np.abs(mu))
        F = (V * w) @ V.conj().T          # d val = Tr(F i dA)
        grad = (D @ (1j * F.T).ravel()).real if P else np.zeros(0)
        return val, grad
    z = np.zeros(P)
    for eps in eps_list:
        if P == 0: break
        r = minimize(f, z, args=(eps,), jac=True, method='L-BFGS-B', options=dict(maxiter=5000, gtol=1e-12, ftol=1e-15))
        z = r.x
    A = Aof(z)
    val = p['II'] - 0.25 * np.sum(np.abs(la.eigvalsh(1j * A)))
    zc = z[:Bc.shape[1]]; zy = z[Bc.shape[1]:]
    return dict(value=float(val), c=cbar + Bc @ zc, y=By @ zy)


def optimize_free_fast(h, n, cfix=None, reflect=True, eps_list=(1e-3, 1e-5, 1e-7, 0.0), maxiter=3000):
    """Sparse-gradient version of optimize_free_lbfgs (same problem, same answer)."""
    from scipy.optimize import minimize
    from common import sum_zero_basis
    p = pauli_coeffs(h)
    m = n - 1
    a = lambda j: 2 * j
    b = lambda j: 2 * j + 1
    bond_entries = [(b, a, 1, -2 * p['XX']), (a, b, 1, 2 * p['YY']), (b, b, 1, -2 * p['XY']), (a, a, 1, 2 * p['YX'])]
    # each entry: (row-site func, col-site func, offset, value), plus on-site Z parts
    def build(c, y):
        A = np.zeros((2 * n, 2 * n))
        j = np.arange(m)
        for (fr, fc, off, val) in bond_entries:
            A[fr(j), fc(j + off)] += c * val
        zs = np.zeros(n); zs[:-1] += c * p['ZI']; zs[1:] += c * p['IZ']; zs += y
        k = np.arange(n)
        A[a(k), b(k)] += -2 * zs
        return A - A.T
    def grads(G):
        # G = d val / dA (antisymmetric part handled by using A - A^T construction): dval = sum G_pq dA_pq
        j = np.arange(m)
        gc = np.zeros(m)
        for (fr, fc, off, val) in bond_entries:
            gc += val * (G[fr(j), fc(j + off)] - G[fc(j + off), fr(j)])
        k = np.arange(n)
        gs = -2 * (G[a(k), b(k)] - G[b(k), a(k)])
        gc += p['ZI'] * gs[:-1] + p['IZ'] * gs[1:]
        return gc, gs
    if cfix is None:
        cbar = np.ones(m) / m; Bc = sum_zero_basis(m, reflect)
    else:
        cbar = np.asarray(cfix, float); Bc = np.zeros((m, 0))
    By = sum_zero_basis(n, reflect)
    nc = Bc.shape[1]
    def unpack(z):
        return cbar + Bc @ z[:nc], By @ z[nc:]
    def f(z, eps):
        c, y = unpack(z)
        A = build(c, y)
        mu, V = la.eigh(1j * A)
        w = mu / np.sqrt(mu ** 2 + eps ** 2) if eps > 0 else np.sign(mu)
        val = np.sum(np.sqrt(mu ** 2 + eps ** 2)) if eps > 0 else np.sum(np.abs(mu))
        F = (V * w) @ V.conj().T          # dval = Tr(F i dA) = sum_pq (i F^T)_pq dA_pq
        G = (1j * F.T).real
        gc, gs = grads(G)
        return val, np.concatenate([Bc.T @ gc, By.T @ gs])
    z = np.zeros(nc + By.shape[1])
    for eps in eps_list:
        r = minimize(f, z, args=(eps,), jac=True, method='L-BFGS-B', options=dict(maxiter=maxiter, gtol=1e-13, ftol=1e-16))
        z = r.x
    c, y = unpack(z)
    val = p['II'] - 0.25 * np.sum(np.abs(la.eigvalsh(1j * build(c, y))))
    return dict(value=float(val), c=c, y=y)

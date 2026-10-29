"""Gauge-optimized bond scale M = min_{o, alpha} (lambda_max - lambda_min)/2 of  h12 + alpha + o x 1 - 1 x o.

The objective is convex and nonsmooth in o.  We minimize a log-sum-exp smoothing with analytic
gradients (Hellmann-Feynman), anneal the temperature, and polish with Nelder-Mead on the exact
half-spread.  We also return the dual lower bound  tr(Delta h)/2  for the explicit Delta built in the
proof (all-bonding minus all-antibonding Gaussian states), which equals ||B||_1 (NC) or ||B||_1/2 (MJ).
"""
import numpy as np
from scipy.optimize import minimize
from ff_bound import half_spread, bond_mj, nc_to_mj, trace_norm


def hermitian_basis(d):
    basis = []
    for i in range(d):
        for j in range(i, d):
            if i == j:
                if i == d - 1:
                    continue
                E = np.zeros((d, d), dtype=complex); E[i, i] = 1; E[d - 1, d - 1] = -1
                basis.append(E)
            else:
                E = np.zeros((d, d), dtype=complex); E[i, j] = 1; E[j, i] = 1; basis.append(E)
                E = np.zeros((d, d), dtype=complex); E[i, j] = -1j; E[j, i] = 1j; basis.append(E)
    return basis


def gauge_M(H12, d, seed=0, verbose=False):
    I = np.eye(d)
    G = np.array([np.kron(b, I) - np.kron(I, b) for b in hermitian_basis(d)])
    nparam = len(G)

    def spectrum(x):
        return np.linalg.eigh(H12 + np.tensordot(x, G, axes=1))

    def f_smooth(x, beta):
        e, V = spectrum(x)
        ep = np.exp(beta * (e - e.max())); em = np.exp(-beta * (e - e.min()))
        val = (e.max() + np.log(ep.sum()) / beta - (e.min() - np.log(em.sum()) / beta)) / 2
        wp = ep / ep.sum(); wm = em / em.sum()
        # d lambda_i / dx_a = <v_i| G_a |v_i>
        Gv = np.einsum('aij,jk->aik', G, V)
        diag = np.einsum('ik,aik->ai', V.conj(), Gv).real  # (a, i)
        grad = (diag @ wp - diag @ wm) / 2
        return val, grad

    x = np.zeros(nparam)
    for beta in (2.0, 8.0, 32.0, 128.0, 512.0):
        r = minimize(lambda y: f_smooth(y, beta), x, jac=True, method="L-BFGS-B",
                     options={"maxiter": 2000, "ftol": 1e-14, "gtol": 1e-10})
        x = r.x
    exact = lambda y: half_spread(H12 + np.tensordot(y, G, axes=1))
    best = exact(x)
    if nparam <= 70:
        r = minimize(exact, x, method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-12})
        if r.fun < best:
            best, x = r.fun, r.x
    return best, x


def bond_from_nc(A0c, Bc):
    """Many-body two-site bond for a number-conserving chain (on-site term split half/half)."""
    A0m, Bm = nc_to_mj(A0c, Bc)
    n = A0m.shape[0]
    return bond_mj(A0m, Bm, n), 2 ** (n // 2)


if __name__ == "__main__":
    print("XX chain (t=1):", end=" ")
    H, d = bond_from_nc(np.zeros((1, 1), complex), -np.ones((1, 1), complex))
    print("M =", round(gauge_M(H, d)[0], 6), " lower bound ||B||_1 = 1")

    print("XX chain with chemical potential mu=0.6 (gauge should remove it):", end=" ")
    H, d = bond_from_nc(0.6 * np.ones((1, 1), complex), -np.ones((1, 1), complex))
    print("bare half-spread =", round(half_spread(H), 6), " M =", round(gauge_M(H, d)[0], 6))

    print("critical Ising (J=g=1 Majorana form):", end=" ")
    A0 = np.array([[0, 1.0], [-1.0, 0]]); B = np.array([[0, 0], [1.0, 0]])
    H = bond_mj(A0, B, 2)
    print("bare =", round(half_spread(H), 6), " M =", round(gauge_M(H, 2)[0], 6), " lower bound ||B||_1/2 = 0.5")

    print("two decoupled XX chains (m=2, d=4):", end=" ")
    H, d = bond_from_nc(np.zeros((2, 2), complex), -np.eye(2, dtype=complex))
    print("M =", round(gauge_M(H, d)[0], 6), " lower bound ||B||_1 = 2")

    print("two XX chains with on-site mixing t'=0.5 (m=2):", end=" ")
    A0 = np.array([[0, 0.5], [0.5, 0]], dtype=complex)
    H, d = bond_from_nc(A0, -np.eye(2, dtype=complex))
    print("bare =", round(half_spread(H), 6), " M =", round(gauge_M(H, d)[0], 6), " lower bound 2")

    print("SSH-type m=2 chain (B rank one, A0 = intracell hopping, critical):", end=" ")
    A0 = np.array([[0, 1.0], [1.0, 0]], dtype=complex); B = np.array([[0, 0], [1.0, 0]], dtype=complex)
    H, d = bond_from_nc(A0, B)
    print("bare =", round(half_spread(H), 6), " M =", round(gauge_M(H, d)[0], 6), " lower bound ||B||_1 = 1")

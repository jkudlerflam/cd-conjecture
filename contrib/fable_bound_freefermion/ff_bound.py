"""
Velocity-weighted Casimir coefficient A versus the bond scale M for
translation-invariant nearest-neighbour quasi-free chains.

Two classes are handled.

Number-conserving (NC): m complex modes per site, d = 2^m,
    H = sum_x [ f_x^dag A0 f_x + f_x^dag B f_{x+1} + h.c. ],
    h(k) = A0 + B e^{ik} + B^dag e^{-ik}  (m x m Hermitian),
    A = (1/2) sum_{Fermi points} |d eps/dk|.

Majorana / Bogoliubov (MJ): n = 2m Majoranas per site,
    H = (i/4) sum_{x,y} gamma_x^T A_{x-y} gamma_y,  A_0 real antisymmetric,
    A_1 = B real, A_{-1} = -B^T,
    h(k) = i (A0 + B e^{ik} - B^T e^{-ik})  (n x n Hermitian),
    A = (1/4) sum_{zero crossings of all n bands} |d eps/dk|.

Fermi points are the unimodular roots of det Q(z), Q(z) = z h(k), found
from the companion linearization of the quadratic matrix polynomial.
Velocities come from Hellmann-Feynman on the kernel of h(k_F)
(degenerate perturbation theory when the kernel has dimension > 1).

M is the gauge-optimized half-spread of the many-body bond on C^d x C^d,
    M = min_{o, alpha} (lambda_max - lambda_min)/2 of  h12 + alpha + o x 1 - 1 x o,
computed by a nonsmooth minimization over the d^2 real parameters of o,
and compared with the proved lower bound ||B||_1 (NC) or ||B||_1 / 2 (MJ).
"""
import numpy as np
from numpy.linalg import eigh, eigvals, svd, norm
from scipy.optimize import minimize


# ----------------------------------------------------------------------
# Fermi points and A
# ----------------------------------------------------------------------
def bloch_nc(A0, B):
    return lambda k: A0 + B * np.exp(1j * k) + B.conj().T * np.exp(-1j * k)


def bloch_mj(A0, B):
    return lambda k: 1j * (A0 + B * np.exp(1j * k) - B.T * np.exp(-1j * k))


def quad_coeffs(kind, A0, B):
    """Q(z) = C0 + C1 z + C2 z^2 with Q(e^{ik}) = e^{ik} h(k)."""
    if kind == "nc":
        return B.conj().T, A0, B
    else:
        return -1j * B.T, 1j * A0, 1j * B


def fermi_points(kind, A0, B, tol=1e-7):
    """Return list of (k_F, kernel basis U) for unimodular roots of det Q."""
    C0, C1, C2 = quad_coeffs(kind, A0, B)
    n = A0.shape[0]
    # generalized companion pencil  (lam*E - F) with E = diag(C2, 1), F = [[-C1, -C0],[1, 0]]
    E = np.block([[C2, np.zeros((n, n))], [np.zeros((n, n)), np.eye(n)]])
    F = np.block([[-C1, -C0], [np.eye(n), np.zeros((n, n))]])
    from scipy.linalg import eig
    w = eig(F, E, right=False)
    w = w[np.isfinite(w)]
    uni = w[np.abs(np.abs(w) - 1.0) < 1e-6]
    ks = np.sort(np.angle(uni))
    # cluster
    pts = []
    for k in ks:
        if pts and abs(np.angle(np.exp(1j * (k - pts[-1])))) < 1e-6:
            continue
        pts.append(k)
    h = bloch_nc(A0, B) if kind == "nc" else bloch_mj(A0, B)
    out = []
    for k in pts:
        e, U = eigh(h(k))
        ker = np.abs(e) < 1e-6 * max(1.0, np.abs(e).max())
        if ker.sum() == 0:
            continue
        out.append((k, U[:, ker]))
    return out


def velocities(kind, A0, B):
    """List of |v| over all crossings (band, k_F), with multiplicity."""
    h = bloch_nc(A0, B) if kind == "nc" else bloch_mj(A0, B)

    def dh(k):
        if kind == "nc":
            return 1j * (B * np.exp(1j * k) - B.conj().T * np.exp(-1j * k))
        return -(B * np.exp(1j * k) + B.T * np.exp(-1j * k))

    vs = []
    for k, U in fermi_points(kind, A0, B):
        V = U.conj().T @ dh(k) @ U
        vs.extend(np.linalg.eigvalsh((V + V.conj().T) / 2).tolist())
    return np.abs(np.array(vs))


def casimir_A(kind, A0, B):
    v = velocities(kind, A0, B)
    return v.sum() / (2.0 if kind == "nc" else 4.0)


# ----------------------------------------------------------------------
# Many-body bond and M
# ----------------------------------------------------------------------
def majoranas(nsites, n):
    """Jordan-Wigner Majoranas gamma_{x,a}, x < nsites, a < n, as 2^(n*nsites/2) matrices."""
    nmodes = nsites * n // 2
    assert nsites * n % 2 == 0
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    I = np.eye(2, dtype=complex)

    def kron_list(ops):
        out = np.array([[1.0 + 0j]])
        for o in ops:
            out = np.kron(out, o)
        return out

    gam = []
    for j in range(nmodes):
        gam.append(kron_list([Z] * j + [X] + [I] * (nmodes - j - 1)))
        gam.append(kron_list([Z] * j + [Y] + [I] * (nmodes - j - 1)))
    # gam[2j], gam[2j+1] are the two Majoranas of complex mode j; mode j = (x, a) ordered site-major
    # reshape to [site][majorana index a]
    out = []
    for x in range(nsites):
        out.append([gam[x * n + a] for a in range(n)])
    return out


def bond_mj(A0, B, n):
    """Two-site bond: (i/2) gamma_1^T B gamma_2 + (i/8)(gamma_1 A0 gamma_1 + gamma_2 A0 gamma_2)."""
    g = majoranas(2, n)
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for a in range(n):
        for b in range(n):
            H += 0.5j * B[a, b] * g[0][a] @ g[1][b]
            H += 0.125j * A0[a, b] * (g[0][a] @ g[0][b] + g[1][a] @ g[1][b])
    return H


def nc_to_mj(A0c, Bc):
    """Complex m-mode (A0, B) -> real Majorana (A0, B) with f = (a + i b)/2 ordering (a1,b1,a2,b2,...)."""
    m = A0c.shape[0]
    n = 2 * m
    # f_a^dag M_ab f_b with M Hermitian  ->  (i/4) gamma^T A gamma + const, A real antisymmetric.
    # f_a = (g_{2a} + i g_{2a+1})/2.  f_a^dag f_b = (1/4)(g_{2a} - i g_{2a+1})(g_{2b} + i g_{2b+1}).
    # do it directly and more carefully
    def quad2(M):
        mm = M.shape[0]
        c = np.zeros((2 * mm, 2 * mm), dtype=complex)
        for a in range(mm):
            for b in range(mm):
                Mab = M[a, b]
                c[2 * a, 2 * b] += 0.25 * Mab
                c[2 * a, 2 * b + 1] += 0.25j * Mab
                c[2 * a + 1, 2 * b] += -0.25j * Mab
                c[2 * a + 1, 2 * b + 1] += 0.25 * Mab
        ca = (c - c.T) / 2
        A = -4j * ca
        assert np.allclose(A.imag, 0)
        return A.real

    # H = sum_x f_x^dag A0c f_x + f_x^dag Bc f_{x+1} + f_{x+1}^dag Bc^dag f_x
    # The on-site part A0c -> A_0;  hopping f_x^dag Bc f_{x+1}  + h.c. -> (i/4)(g_x^T B g_{x+1} - g_{x+1}^T B^T g_x)
    # Build via the 2-site quadratic form on 2m complex modes, then read off blocks.
    Mbig = np.zeros((2 * m, 2 * m), dtype=complex)
    Mbig[:m, m:] = Bc
    Mbig[m:, :m] = Bc.conj().T
    Abig = quad2(Mbig)  # (2n x 2n) real antisymmetric
    Bmj = Abig[:n, n:]
    assert np.allclose(Abig[n:, :n], -Bmj.T)
    A0mj = quad2(A0c)
    return A0mj, Bmj


def half_spread(H):
    e = np.linalg.eigvalsh(H)
    return (e[-1] - e[0]) / 2


def gauge_opt_M(H12, d, restarts=3, seed=0):
    """min over Hermitian o of half-spread of H12 + o x 1 - 1 x o."""
    rng = np.random.default_rng(seed)
    I = np.eye(d)
    # Hermitian basis of traceless o: d^2 - 1 real params (trace drops out)
    basis = []
    for i in range(d):
        for j in range(i, d):
            E = np.zeros((d, d), dtype=complex)
            if i == j:
                if i == d - 1:
                    continue
                E[i, i] = 1
                E[d - 1, d - 1] = -1
                basis.append(E)
            else:
                E[i, j] = 1
                E[j, i] = 1
                basis.append(E.copy())
                E2 = np.zeros((d, d), dtype=complex)
                E2[i, j] = -1j
                E2[j, i] = 1j
                basis.append(E2)
    G = np.array([np.kron(b, I) - np.kron(I, b) for b in basis])

    def f(x):
        return half_spread(H12 + np.tensordot(x, G, axes=1))

    best = half_spread(H12)
    x0s = [np.zeros(len(basis))] + [0.3 * rng.standard_normal(len(basis)) for _ in range(restarts - 1)]
    for x0 in x0s:
        r = minimize(f, x0, method="Nelder-Mead", options={"maxiter": 20000, "xatol": 1e-9, "fatol": 1e-11})
        r = minimize(f, r.x, method="Powell", options={"maxiter": 20000, "xtol": 1e-9, "ftol": 1e-12})
        best = min(best, r.fun)
    return best


# ----------------------------------------------------------------------
# Random sampling
# ----------------------------------------------------------------------
def random_nc(m, rng, onsite=True):
    B = rng.standard_normal((m, m)) + 1j * rng.standard_normal((m, m))
    A0 = np.zeros((m, m), dtype=complex)
    if onsite:
        X = rng.standard_normal((m, m)) + 1j * rng.standard_normal((m, m))
        A0 = (X + X.conj().T) / 2
    return A0, B


def random_mj(n, rng, onsite=True):
    B = rng.standard_normal((n, n))
    A0 = np.zeros((n, n))
    if onsite:
        X = rng.standard_normal((n, n))
        A0 = (X - X.T) / 2
    return A0, B


def trace_norm(B):
    return svd(B, compute_uv=False).sum()


if __name__ == "__main__":
    import sys
    rng = np.random.default_rng(1)
    print("=== sanity: XX chain at half filling, m=1 ===")
    A0 = np.zeros((1, 1), dtype=complex); B = -np.ones((1, 1), dtype=complex)
    print("A (nc) =", casimir_A("nc", A0, B), " ||B||_1 =", trace_norm(B))
    A0m, Bm = nc_to_mj(A0, B)
    print("Majorana B =\n", Bm, "\nA (mj) =", casimir_A("mj", A0m, Bm), " ||B_mj||_1/2 =", trace_norm(Bm) / 2)
    H12 = bond_mj(A0m, Bm, 2)
    print("bond spectrum", np.round(np.linalg.eigvalsh(H12), 6), " gauge-opt M =", gauge_opt_M(H12, 2))

    print("\n=== sanity: critical Ising, n=2 Majoranas ===")
    A0 = np.array([[0, 1.0], [-1.0, 0]]); B = np.array([[0, 0], [1.0, 0]])
    print("A (mj) =", casimir_A("mj", A0, B), " ||B||_1/2 =", trace_norm(B) / 2)
    H12 = bond_mj(A0, B, 2)
    print("bare half-spread", half_spread(H12), " gauge-opt M =", gauge_opt_M(H12, 2))

    print("\n=== random number-conserving chains: A / ||B||_1 ===")
    for m in (1, 2, 3, 4):
        ratios = []
        for _ in range(2000):
            A0, B = random_nc(m, rng)
            ratios.append(casimir_A("nc", A0, B) / trace_norm(B))
        print(f"m={m}: max A/||B||_1 = {max(ratios):.4f}  (bound claims 2)")

    print("\n=== random Majorana chains: A / (||B||_1/2) ===")
    for n in (1, 2, 3, 4, 6):
        ratios = []
        for _ in range(2000):
            A0, B = random_mj(n, rng)
            ratios.append(casimir_A("mj", A0, B) / (trace_norm(B) / 2))
        print(f"n={n}: max A/(||B||_1/2) = {max(ratios):.4f}  (bound claims 2)")

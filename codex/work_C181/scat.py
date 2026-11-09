"""C175: exact lattice two-body threshold scattering for the Fibonacci-cusp model (C173 PROOF eq. 3.2).

Model (units |t2| = 1): hard-core particles, no two adjacent, H = sum_j [w c+_{j+1} c_j + h.c. - (c+_{j+2} c_j + h.c.)
 - 2 n_j n_{j+2} - (2/3) n_j n_{j+3}].  Single-particle hop amplitudes a_d (amplitude to move by +d):
a_1 = w, a_-1 = conj(w), a_{+-2} = t2 = -1.  eps(k) = sum_d a_d e^{-ikd} = 2|w|cos(k - phi) - 2 cos 2k.

Two particles x1 < x2, r = x2 - x1 >= 2 (blockade).  psi(x1,x2) = e^{iK x1} f(r):
  (h f)(r) = sum_d a_d [ e^{-iKd} f(r+d) + f(r-d) ] + V(r) f(r),  f(r) = 0 for r < 2.
A plane wave f = e^{ipr} has right momentum p, left momentum K - p, energy eps(K-p) + eps(p).

Threshold data.  Free modes of the 4th-order recurrence at energy E are z^r with z roots of a quartic.
At K = k_a + k_b (= pi at phi = pi/2) and E = 2 eps_min the roots are double at e^{ik_a}, e^{ik_b}:
  f(r) = e^{ik_b r}(alpha_ab + beta_ab r) + e^{ik_a r}(alpha_ba + beta_ba r)   ("ab" = left a, right b)
Two boundary equations (r = 2, 3) leave a 2-dim space; alpha = -A beta defines the 2x2 matrix A.
At K = 2 k_a: one double root (open channel) and two simple roots (one decaying, one growing); keep the decaying one.
"""
import numpy as np

def hops(w, t2=-1.0):
    return {1: w, -1: np.conj(w), 2: t2, -2: t2}

def eps(k, w, t2=-1.0):
    return sum(a*np.exp(-1j*k*d) for d, a in hops(w, t2).items()).real

def valleys(w, t2=-1.0, n=200001):
    ks = np.linspace(-np.pi, np.pi, n, endpoint=False)
    e = eps(ks, w, t2)
    # local minima on the circle
    idx = np.nonzero((e < np.roll(e, 1)) & (e < np.roll(e, -1)))[0]
    out = []
    for i in idx:
        k = ks[i]
        for _ in range(60):  # Newton polish on eps'
            h = 1e-5
            d1 = (eps(k+h, w, t2)-eps(k-h, w, t2))/(2*h)
            d2 = (eps(k+h, w, t2)-2*eps(k, w, t2)+eps(k-h, w, t2))/h**2
            k -= d1/d2
        out.append((k, eps(k, w, t2)))
    out.sort(key=lambda x: x[1])
    return out

def rel_coeffs(K, w, t2=-1.0):
    """Coefficient c[s] multiplying f(r+s) in (h f)(r), s = -2..2 (without V and E)."""
    a = hops(w, t2)
    c = {s: 0j for s in (-2, -1, 0, 1, 2)}
    for d, ad in a.items():
        c[d] += ad*np.exp(-1j*K*d)   # left particle moves by -d -> r+d
        c[-d] += ad                  # right particle moves by -d -> r-d
    return c

def free_roots(K, E, w, t2=-1.0):
    c = rel_coeffs(K, w, t2)
    # sum_s c[s] z^s - E = 0  ->  multiply by z^2
    poly = [c[2], c[1], c[0]-E, c[-1], c[-2]]
    return np.roots(poly)

def Vpot(r, V2=-2.0, V3=-2/3, V1hard=True):
    return {2: V2, 3: V3}.get(r, 0.0)

def boundary_rows(modes, K, E, w, t2=-1.0, V2=-2.0, V3=-2/3, rmin=2, nnn=True):
    """Rows of the homogeneous boundary equations at r = rmin, rmin+1 for coefficients of the free modes.
    modes: list of functions r -> value.  True eq at r: free eq minus terms with f(r') for r' < rmin, plus V(r) f(r)."""
    c = rel_coeffs(K, w, t2)
    rows = []
    for r in (rmin, rmin+1):
        row = []
        for m in modes:
            val = 0j
            for s in (-2, -1):
                if r+s < rmin:
                    val -= c[s]*m(r+s)
            val += Vpot(r, V2, V3)*m(r)
            row.append(val)
        rows.append(row)
    return np.array(rows)

def A_matrix(w, t2=-1.0, V2=-2.0, V3=-2/3, E=None, check=True):
    """2x2 threshold scattering-length matrix for the unlike pair at K = k_a + k_b. Channel order (ab, ba)."""
    (ka, ea), (kb, eb) = valleys(w, t2)[:2]
    K = ka + kb
    if E is None: E = ea + eb
    modes = [lambda r, k=kb: np.exp(1j*k*r), lambda r, k=ka: np.exp(1j*k*r),
             lambda r, k=kb: r*np.exp(1j*k*r), lambda r, k=ka: r*np.exp(1j*k*r)]
    M = boundary_rows(modes, K, E, w, t2, V2, V3)
    Ma, Mb = M[:, :2], M[:, 2:]
    A = np.linalg.solve(Ma, Mb)   # Ma alpha + Mb beta = 0  ->  alpha = -A beta
    if check:
        # free equation check for the modes (double roots)
        c = rel_coeffs(K, w, t2)
        for m in modes:
            res = sum(c[s]*m(10+s) for s in c) - E*m(10)
            assert abs(res) < 1e-6, res
    return A, ka, kb, K

def a_like(w, valley=0, t2=-1.0, V2=-2.0, V3=-2/3):
    """Threshold scattering length of a like pair in valley 0 (lower or a) or 1."""
    k, e = valleys(w, t2)[valley]
    K = 2*k; E = 2*e
    z = free_roots(K, E, w, t2)
    zk = np.exp(1j*k)
    others = sorted(z, key=lambda x: abs(x-zk))[2:]
    # double root sanity
    dec = [x for x in others if 1e-12 < abs(x) < 1-1e-9]
    opn = [x for x in others if abs(abs(x)-1) <= 1e-9]
    assert not opn, ("extra open channel", others)
    modes = [lambda r: zk**r, lambda r: r*zk**r] + [lambda r, x=x: x**r for x in dec]
    M = boundary_rows(modes, K, E, w, t2, V2, V3)
    M = M[np.abs(M).max(axis=1) > 1e-13]
    # unknowns (alpha, beta, gamma...): solve for alpha, gammas in terms of beta = 1
    cols = [0] + list(range(2, M.shape[1]))
    sol = np.linalg.solve(M[:, cols], -M[:, 1])
    alpha = sol[0]
    return -alpha, k  # f ~ e^{ikr}(r - a) * beta  ->  a = -alpha/beta

def bound_scan(w, t2=-1.0, V2=-2.0, V3=-2/3, Lr=600, nK=97):
    """Lowest two-body level minus continuum bottom, over total momentum K (open relative box of Lr)."""
    import scipy.linalg as sl
    ks = np.linspace(-np.pi, np.pi, 4001)
    worst = (np.inf, None)
    res = []
    for K in np.linspace(-np.pi, np.pi, nK):
        c = rel_coeffs(K, w, t2)
        n = Lr
        # band matrix in r = 2..Lr+1
        H = np.zeros((n, n), complex)
        for i in range(n):
            r = i+2
            H[i, i] = c[0] + Vpot(r, V2, V3)
            for s in (1, 2):
                if i+s < n:
                    H[i, i+s] = c[s]
                    H[i+s, i] = c[-s]
        e0 = sl.eigh(H, eigvals_only=True, subset_by_index=[0, 0])[0]
        cont = (eps(ks, w, t2)+eps(K-ks, w, t2)).min()
        res.append((K, e0-cont))
        if e0-cont < worst[0]: worst = (e0-cont, K)
    return worst, res

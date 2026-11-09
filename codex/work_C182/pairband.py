"""C182: exact two-body bands of the Fibonacci-cusp model (C173 PROOF eq. 3.2), c frame, units |t2| = 1.

H = sum_j [ w c+_{j+1} c_j + h.c. - (c+_{j+2} c_j + h.c.) - 2 n_j n_{j+2} - (2/3) n_j n_{j+3} ], blockade.
Single-particle hop a_d (amplitude to move by +d): a_1 = w, a_-1 = conj(w), a_{+-2} = t2 = -1.
eps(k) = 2|w| cos(k - phi) - 2 cos 2k.

Two particles x1 < x2, r = x2 - x1 >= 2, psi = e^{iK x1} f(r):
  (h_K f)(r) = sum_s c_s f(r+s) + V(r) f(r),  c_d = a_d e^{-iKd} + a_{-d}  (same convention as C175 scat.py).
The half-line operator truncated to 2 <= r <= R+1 is a Dirichlet restriction, so its eigenvalues are upper
bounds (min-max). E_0(K; R) < E_cont(K) therefore certifies a bound state at K.
"""
import numpy as np, scipy.linalg as sl, sys, json

def hops(w, t2=-1.0):
    return {1: w, -1: np.conj(w), 2: t2, -2: t2}

def eps(k, w, t2=-1.0):
    return (2*np.real(np.conj(w)*np.exp(1j*k)) + 2*t2*np.cos(2*k)) if True else None

def eps_check(k, w, t2=-1.0):
    return sum(a*np.exp(-1j*k*d) for d, a in hops(w, t2).items()).real

def valleys(w, t2=-1.0, n=20001):
    ks = np.linspace(-np.pi, np.pi, n, endpoint=False)
    e = eps(ks, w, t2)
    idx = np.nonzero((e < np.roll(e, 1)) & (e < np.roll(e, -1)))[0]
    out = []
    for i in idx:
        k = ks[i]
        for _ in range(50):
            h = 1e-5
            d1 = (eps(k+h, w, t2)-eps(k-h, w, t2))/(2*h)
            d2 = (eps(k+h, w, t2)-2*eps(k, w, t2)+eps(k-h, w, t2))/h**2
            k -= d1/d2
        out.append((float(np.angle(np.exp(1j*k))), float(eps(k, w, t2))))
    out.sort(key=lambda x: x[1])
    return out

_KS = np.linspace(-np.pi, np.pi, 8192, endpoint=False)
def econt(K, w, t2=-1.0):
    """Bottom of the two-particle scattering continuum at total momentum K (grid + golden refinement)."""
    e = eps(_KS, w, t2) + eps(K-_KS, w, t2)
    i = int(np.argmin(e)); a, b = _KS[i]-1e-3, _KS[i]+1e-3
    f = lambda p: eps(p, w, t2)+eps(K-p, w, t2)
    g = (np.sqrt(5)-1)/2
    for _ in range(60):
        c, d = b-g*(b-a), a+g*(b-a)
        if f(c) < f(d): b = d
        else: a = c
    return min(e[i], f(0.5*(a+b)))

def rel_band(K, w, R, t2=-1.0, V2=-2.0, V3=-2/3):
    a = hops(w, t2)
    c = {s: 0j for s in (-2, -1, 0, 1, 2)}
    for d, ad in a.items():
        c[d] += ad*np.exp(-1j*K*d)
        c[-d] += ad
    n = R
    ab = np.zeros((3, n), complex)          # upper form, u = 2: ab[2 + i - j, j] = A[i, j]
    diag = np.full(n, c[0].real)
    diag[0] += V2; diag[1] += V3            # r = 2, 3
    ab[2, :] = diag
    ab[1, 1:] = c[1]                        # A[i, i+1] = c_1  (f(r+1) coefficient in row r)
    ab[0, 2:] = c[2]
    return ab

def lowest(K, w, R, k=2, **kw):
    ab = rel_band(K, w, R, **kw)
    return sl.eig_banded(ab, lower=False, eigvals_only=True, select='i', select_range=(0, k-1))

def scan(wa, phi, R=3000, nK=361, **kw):
    w = wa*np.exp(1j*phi)
    Ks = np.linspace(-np.pi, np.pi, nK, endpoint=False)
    E0 = np.array([lowest(K, w, R, k=1, **kw)[0] for K in Ks])
    Ec = np.array([econt(K, w) for K in Ks])
    return Ks, E0, Ec

def Eb_exact(u):
    """ab pair binding at phi = pi/2, collaborator's closed form (contrib/cusp_scattering_20261109/DERIVATION.md eqs 1, 2)."""
    from scipy.optimize import brentq
    p = brentq(lambda p: (1-p)*(3-p)*(1+p)**2 - u*u*p**3, 1e-12, 1-1e-12)
    Ep = -2*(p+1/p+u*u*p/(1+p)**2)
    return (-4-u*u/2) - Ep

def refine_min(f, K0, h0):
    a, b = K0-h0, K0+h0
    g = (np.sqrt(5)-1)/2
    for _ in range(40):
        c, d = b-g*(b-a), a+g*(b-a)
        if f(c) < f(d): b = d
        else: a = c
    K = 0.5*(a+b); return K, f(K)

def analyse(wa, phi, R=3000, nK=361, **kw):
    w = wa*np.exp(1j*phi)
    Ks, E0, Ec = scan(wa, phi, R, nK, **kw)
    v = valleys(w)
    ea = v[0][1]
    dK = Ks[1]-Ks[0]
    # local minima of the lower edge min(E0, Ec) (E0 <= Ec up to box error anyway)
    low = np.minimum(E0, Ec)
    mins = []
    idx = np.nonzero((low <= np.roll(low, 1)) & (low <= np.roll(low, -1)))[0]
    for i in idx:
        f = lambda K: min(lowest(K, w, R, k=1, **kw)[0], econt(K, w))
        K, e = refine_min(f, Ks[i], dK)
        ec = econt(K, w)
        e2 = lowest(K, w, 2*R, k=1, **kw)[0]
        # curvature of the bound band
        h = 2e-3
        curv = (f(K+h)-2*f(K)+f(K-h))/h**2
        mins.append(dict(K=float(np.angle(np.exp(1j*K))), E=float(e), E_2R=float(e2), depth=float(ec-e),
                         E_minus_2ea=float(e-2*ea), curv=float(curv)))
    mins.sort(key=lambda m: m['E'])
    return dict(wa=wa, phi=phi, R=R, valleys=v, delta=(v[1][1]-v[0][1]) if len(v) > 1 else None,
                mins=mins, Ks=Ks.tolist(), bound=(Ec-E0).tolist())

if __name__ == "__main__":
    # self-checks
    w = 1j*1.0
    ks = np.linspace(-3, 3, 7)
    assert np.allclose(eps(ks, w), eps_check(ks, w))
    # collaborator's exact value at u = 1: E_pair = -4.519526534666399 (real gauge K = 0 <-> c frame K = pi)
    e = lowest(np.pi, w, 4000, k=1)[0]
    print("check |w|=1 phi=pi/2 K=pi: E0 =", repr(e), " exact -4.519526534666399", flush=True)
    e = lowest(np.pi, 2j, 4000, k=1)[0]
    print("check |w|=2: E0 =", repr(e), " exact -6.122547914236252", flush=True)

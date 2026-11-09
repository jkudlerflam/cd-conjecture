"""Check of PROOF §5 (C169): free (JW-quadratic) nearest-neighbour qubit chains of any period p.

Majoranas a_j, b_j on site j.  Allowed terms (Pauli -> Majorana bilinear):
  Z_j -> a_j b_j,  X_j X_{j+1} -> b_j a_{j+1},  Y_j Y_{j+1} -> a_j b_{j+1},
  X_j Y_{j+1} -> b_j b_{j+1},  Y_j X_{j+1} -> a_j a_{j+1}.
Bloch matrix h(z) of size 2p; the claim is that z^2 det h(z) is a polynomial of
degree <= 4 (<= 2 if one bond has rank-1 in-plane block), so at most 4 zero
crossings on |z| = 1 and c <= 1.  We verify the degree by exact sampling on a
circle of radius r != 1 and count unit-circle roots for random chains, and
compare with range-2 chains (degree <= 8, c <= 2).
"""
import numpy as np
rng = np.random.default_rng(1)

def bloch(D, T, z):
    """D: list of p on-site couplings, T: list of p 2x2 inter-site blocks
    (rows a_j,b_j ; cols a_{j+1},b_{j+1}).  Returns A(z) (antisymmetric part
    convention A(z)^T = -A(1/z))."""
    p = len(D)
    A = np.zeros((2 * p, 2 * p), complex)
    for j in range(p):
        A[2*j, 2*j+1] += D[j]; A[2*j+1, 2*j] -= D[j]
        k = (j + 1) % p
        ph = z if k == 0 else 1.0
        A[2*j:2*j+2, 2*k:2*k+2] += T[j] * ph
        A[2*k:2*k+2, 2*j:2*j+2] -= T[j].T / ph
    return A

def poly_degree(f, nmax=12, r=1.0):
    # sample g(z) = z^nmax/2 * f(z) on a circle and read Laurent coefficients
    N = 64
    zs = r * np.exp(2j * np.pi * np.arange(N) / N)
    vals = np.array([f(z) for z in zs])
    c = np.fft.fft(vals) / N            # c[m] ~ coefficient of z^m * r^m (m mod N)
    m = np.arange(N); m[m > N // 2] -= N
    coef = c / r**m
    nz = [mm for mm, cc in zip(m, coef) if abs(cc) > 1e-9 * max(1, abs(coef).max())]
    return min(nz), max(nz)

def unit_roots(D, T):
    f = lambda z: np.linalg.det(bloch(D, T, z))
    lo, hi = poly_degree(f)
    # polynomial z^{-lo} f(z) of degree hi-lo; get coefficients exactly by fft
    N = 64; r = 1.0
    zs = np.exp(2j*np.pi*np.arange(N)/N)
    c = np.fft.fft([f(z) for z in zs]) / N
    m = np.arange(N); m[m > N//2] -= N
    coeffs = {mm: cc for mm, cc in zip(m, c)}
    P = [coeffs[k] for k in range(hi, lo - 1, -1)]
    roots = np.roots(P)
    return lo, hi, int(np.sum(np.abs(np.abs(roots) - 1) < 1e-6))

if __name__ == '__main__':
    worst = 0
    for p in (1, 2, 3, 4):
        for trial in range(200):
            D = list(rng.normal(size=p))
            T = [rng.normal(size=(2, 2)) for _ in range(p)]
            lo, hi = poly_degree(lambda z: np.linalg.det(bloch(D, T, z)))
            worst = max(worst, hi - lo)
        print(f'p={p}: max Laurent span of det h(z) over 200 random chains = {worst}')
    # rank-1 bond somewhere in the cell
    for p in (2, 3):
        D = list(rng.normal(size=p)); T = [rng.normal(size=(2, 2)) for _ in range(p)]
        u, v = rng.normal(size=2), rng.normal(size=2); T[0] = np.outer(u, v)
        lo, hi = poly_degree(lambda z: np.linalg.det(bloch(D, T, z)))
        print(f'p={p}, one rank-1 bond: span {hi-lo}')
    # examples: XX chain (c=1), TFI critical (c=1/2), Kitaev chain = XZ point (flat band)
    XX = np.array([[0, 1.], [-1., 0]])  # XX -> -i b_j a_{j+1}, YY -> +i a_j b_{j+1}: opposite signs
    print('XX chain mu=0.3:', unit_roots([0.3], [XX]), '(4 unit roots -> c = 1)')
    TFI = np.array([[0, 0], [1., 0]])  # b_j a_{j+1}
    print('TFI critical:', unit_roots([1.0], [TFI]), '(2 roots, double at z=1 -> c = 1/2)')
    # search: maximize number of unit roots over random period-3 chains
    best = 0
    for trial in range(3000):
        p = rng.integers(1, 5)
        D = list(rng.normal(size=p)); T = [rng.normal(size=(2, 2)) for _ in range(p)]
        try:
            best = max(best, unit_roots(D, T)[2])
        except Exception:
            pass
    print('max number of unit-circle roots over 3000 random chains, p <= 4:', best)

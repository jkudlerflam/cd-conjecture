"""C185 R4: first-order dynamics at the Fibonacci point.
Project every translation-invariant bond Pauli sum onto the blockaded space P (no two adjacent particles,
particle = spin down, n = (1-Z)/2) on a ring, and decompose into the first-order operators
  T  = sum_j P(s+_j s-_{j+1} + h.c.)P   (real NN hop)      Ti = sum_j P(i s+_j s-_{j+1} + h.c.)P
  OX = sum_j P X_j P (PXP)                                  OY = sum_j P Y_j P
  N  = sum_j n_j ,  1.
Residual norm printed: zero residual means the projected term lies in this span.
"""
import numpy as np, itertools
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex); sp = np.array([[0, 1], [0, 0]], complex)  # s+ raises: |dn> -> |up>  (|up> = index 0)
sm = sp.T.copy()
paul = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}

def op(L, ops):
    M = np.array([[1.0+0j]])
    for j in range(L):
        M = np.kron(M, ops.get(j, I2))
    return M

def ring_sum(L, A, B):
    return sum(op(L, {j: A, (j+1) % L: B}) for j in range(L)) if B is not None else sum(op(L, {j: A}) for j in range(L))

def Pidx(L):
    out = []
    for s in range(2**L):
        occ = [(s >> (L-1-j)) & 1 for j in range(L)]
        if not any(occ[j] and occ[(j+1) % L] for j in range(L)):
            out.append(s)
    return out

if __name__ == "__main__":
    L = 10
    P = Pidx(L)
    pr = lambda M: M[np.ix_(P, P)]
    # particle = down = basis index 1. s- (spin lowering) creates a particle.
    cre, ann = sm, sp
    T = pr(ring_sum(L, cre, ann) + ring_sum(L, ann, cre))
    Ti = pr(1j*ring_sum(L, cre, ann) - 1j*ring_sum(L, ann, cre))
    OX = pr(ring_sum(L, X, None)); OY = pr(ring_sum(L, Y, None))
    N = pr(ring_sum(L, (I2 - Z)/2, None)); One = np.eye(len(P))
    basis = {'T': T, 'Ti': Ti, 'OX': OX, 'OY': OY, 'N': N, '1': One}
    names = list(basis)
    Bm = np.array([basis[k].ravel() for k in names]).T
    print(f"ring L = {L}, dim P = {len(P)}. Projected bond sum = sum of coefficients x basis; resid = residual norm / norm")
    rows = []
    for a, b in itertools.product('IXYZ', repeat=2):
        if a == 'I' and b == 'I':
            continue
        M = pr(ring_sum(L, paul[a], paul[b]))
        c, *_ = np.linalg.lstsq(Bm, M.ravel(), rcond=None)
        res = np.linalg.norm(Bm @ c - M.ravel()) / max(np.linalg.norm(M), 1e-12)
        cs = '  '.join(f"{k}:{c[i].real:+.3f}" + (f"{c[i].imag:+.3f}i" if abs(c[i].imag) > 1e-9 else '') for i, k in enumerate(names) if abs(c[i]) > 1e-9)
        print(f"  {a}{b}: {cs:60s} resid {res:.1e}")

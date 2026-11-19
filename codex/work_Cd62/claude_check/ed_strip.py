"""Independent ED check of the strip bound L(E_2L - 2E_L)/M <= 2 pi sin(eta) + pi^2/4.

Random U(1) qubit bonds h = Jp(XX+YY)/4 + D(XY-YX)/4 + Jz ZZ/4 - (B/2)(Z1+1Z)/2 with
sqrt(Jp^2+D^2) = 1 and Jz = -cos(eta), 0 < eta <= pi/4.  Since M >= J/2 we test the
stronger quantity L(E_2L - 2E_L)/(J/2).  Also report the gauge-optimized M numerically.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.optimize import minimize

sx = np.array([[0, 1], [1, 0]]) / 2
sy = np.array([[0, -1j], [1j, 0]]) / 2
sz = np.array([[1, 0], [0, -1]]) / 2
I2 = np.eye(2)


def bond(Jp, D, Jz, B):
    return (Jp * (np.kron(sx, sx) + np.kron(sy, sy))
            + D * (np.kron(sx, sy) - np.kron(sy, sx))
            + Jz * np.kron(sz, sz)
            - (B / 2) * (np.kron(sz, I2) + np.kron(I2, sz)))


def ring_E0(h, L):
    """Ground energy of sum_x h_{x,x+1} on an L-site ring (sparse, magnetization blocks via full space)."""
    dim = 2 ** L
    H = sp.csr_matrix((dim, dim), dtype=complex)
    hs = sp.csr_matrix(h)
    for x in range(L - 1):
        H = H + sp.kron(sp.kron(sp.identity(2 ** x), hs), sp.identity(2 ** (L - x - 2)))
    # bond (L-1, 0): permute h so that it acts on last and first site
    h4 = h.reshape(2, 2, 2, 2)  # out1 out2 in1 in2 with site1=L-1, site2=0
    # operator on site0 (first factor) and site L-1 (last factor): swap roles
    hsw = h4.transpose(1, 0, 3, 2).reshape(4, 4)  # now first factor = site 0, second = site L-1
    # build as sum of products via operator Schmidt decomposition
    M = h4.transpose(1, 3, 0, 2).reshape(4, 4)  # (out0,in0),(outL-1,inL-1)
    U, s, Vh = np.linalg.svd(M)
    for k in range(4):
        if s[k] < 1e-14:
            continue
        A = (U[:, k] * s[k]).reshape(2, 2)  # site 0
        Bm = Vh[k, :].reshape(2, 2)         # site L-1
        H = H + sp.kron(sp.kron(sp.csr_matrix(A), sp.identity(2 ** (L - 2))), sp.csr_matrix(Bm))
    if dim <= 512:
        return np.linalg.eigvalsh(H.toarray()).min()
    return sla.eigsh(H, k=1, which='SA', tol=1e-12)[0][0]


def gauge_M(h):
    def spread(p):
        o = np.array([[p[0], p[1] + 1j * p[2]], [p[1] - 1j * p[2], p[3]]])
        hh = h + np.kron(o, I2) - np.kron(I2, o)
        e = np.linalg.eigvalsh(hh)
        return (e[-1] - e[0]) / 2
    best = min(minimize(spread, np.random.randn(4) * 0.3, method='Nelder-Mead',
                        options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 4000}).fun for _ in range(6))
    return best


rng = np.random.default_rng(7)
worst = 0.0
rows = []
for trial in range(24):
    eta = rng.uniform(0.02, np.pi / 4)
    phi = rng.uniform(0, 2 * np.pi)
    Jp, D = np.cos(phi), np.sin(phi)
    B = rng.uniform(-1.2, 1.2) * (1 - np.cos(eta)) * 2  # fields of the scale that keep the sector nontrivial
    h = bond(Jp, D, -np.cos(eta), B)
    Mg = gauge_M(h)
    bound = 2 * np.pi * np.sin(eta) + np.pi ** 2 / 4
    for L in range(2, 8):
        EL = ring_E0(h, L)
        E2L = ring_E0(h, 2 * L)
        q = L * (E2L - 2 * EL) / 0.5
        worst = max(worst, q / bound)
        rows.append((trial, L, eta, B, q, bound, L * (E2L - 2 * EL) / Mg, Mg))
for r in rows[:: max(1, len(rows) // 24)]:
    print("trial %2d L=%d eta=%.3f B=%+.3f  L*def/(J/2)=%+.4f  bound=%.4f  L*def/M=%+.4f  M=%.4f" % r)
print("checked", len(rows), "cases; largest L*def/(J/2) / bound =", round(worst, 4))
print("smallest gauge M over bonds:", round(min(r[7] for r in rows), 4), "(proof uses M >= J/2 = 0.5)")

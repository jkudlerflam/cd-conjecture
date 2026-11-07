"""Exact free-fermion ground energy and entanglement for open qubit chains in the free XZ neighbourhood.
JW along y: beta_j = c_{2j-1} = S_j X_j, alpha_j = c_{2j} = S_j Z_j, S_j = prod_{k<j} Y_k. Dictionary (PROOF 1.1):
  XZ = i a_j a_{j+1}, ZX = -i b_j b_{j+1}, XX = i a_j b_{j+1}, ZZ = -i b_j a_{j+1}, Y_j = i b_j a_j
(a = alpha, b = beta). Coefficients: dict with keys xz, zx, xx, zz, By (bond coefficient; field per site)."""
import numpy as np

def A_matrix(L, cf, open_chain=True):
    n = 2 * L
    A = np.zeros((n, n))
    al = lambda j: 2 * j + 1   # alpha_j (j = 0..L-1) -> c index 2j+1 (0-based c_{2j})
    be = lambda j: 2 * j       # beta_j -> c_{2j-1}
    def add(m, k, c):          # term i c gamma_m gamma_k
        A[m, k] += c; A[k, m] -= c
    nb = L - 1 if open_chain else L
    for j in range(nb):
        k = (j + 1) % L
        add(al(j), al(k), cf.get("xz", 0.0))
        add(be(j), be(k), -cf.get("zx", 0.0))
        add(al(j), be(k), cf.get("xx", 0.0))
        add(be(j), al(k), -cf.get("zz", 0.0))
    for j in range(L):
        wgt = 0.5 if (open_chain and j in (0, L - 1)) else 1.0   # bond fields: half weight at open ends
        add(be(j), al(j), wgt * cf.get("By", 0.0))
    return A

def ground(L, cf):
    A = A_matrix(L, cf)
    h = 2j * A
    w, U = np.linalg.eigh(h)
    E0 = -0.5 * w[w > 0].sum()
    P = U[:, w > 1e-12] @ U[:, w > 1e-12].conj().T
    G = 2 * P - np.eye(2 * L)
    S = []
    for l in range(1, L):
        nu = np.linalg.eigvalsh(G[:2 * l, :2 * l])
        nu = np.clip(nu[nu > 1e-12], 0, 1 - 1e-15)
        p = (1 + nu) / 2
        S.append(float(np.sum(-p * np.log(p) - (1 - p) * np.log(1 - p))))
    return E0, np.array(S)

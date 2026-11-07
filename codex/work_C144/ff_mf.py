"""Reference entanglement of the mean-field state on the tuned line: A dimerized (hoppings 1 +- d_A on XZ bonds),
B staggered (ZX coefficients t (-1)^j + u). Open chain, exact free fermions. Fits with cfit.fits."""
import numpy as np
from cfit import fits
def run(L, dA, t, u):
    n = 2 * L; A = np.zeros((n, n))
    al = lambda j: 2 * j + 1; be = lambda j: 2 * j
    def add(m, k, c): A[m, k] += c; A[k, m] -= c
    for j in range(L - 1):
        add(al(j), al(j + 1), 1 + dA * (-1) ** j)
        add(be(j), be(j + 1), -(u + t * (-1) ** j))
    h = 2j * A; w, U = np.linalg.eigh(h)
    P = U[:, w > 1e-12] @ U[:, w > 1e-12].conj().T
    G = 2 * P - np.eye(n); S = []
    for l in range(1, L):
        nu = np.linalg.eigvalsh(G[:2 * l, :2 * l]); nu = np.clip(nu[nu > 1e-12], 0, 1 - 1e-15); p = (1 + nu) / 2
        S.append(float(np.sum(-p * np.log(p) - (1 - p) * np.log(1 - p))))
    return np.array(S)
for L in (64, 128, 256):
    for u in (0.0, 0.005):
        S = run(L, 0.19, 0.08, u)
        print(f"L={L} dA=0.19 t_s=0.08 u={u}: S_mid={S[L//2-1]:.3f} c(w)={fits(S, L)}")

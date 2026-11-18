"""(B2) check: exchange class at V=0, L(L e0 - E_L) <= 2 t L (csc(pi/L) - L/pi) <= (4-8/pi) t, with random W and mu.
Control outside the class: SU(3) swap chain (three mutually exchanging species) exceeds the (B2) constant.
Claude Fable 5.2, work_C201."""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import unitary_group
import scipy.sparse as sp, scipy.sparse.linalg as spl
from exchange_class_check import bond, chain_H
rng = np.random.default_rng(11)
worst = 0
for (d0, d1) in ((1, 1), (2, 1), (1, 2), (2, 2)):
    d = d0+d1
    for mu in (0.0, 0.3, -0.7, 1.2, 1.9, 2.0, -2.0, 2.5):
        W = unitary_group.rvs(d0*d1, random_state=rng) if d0*d1 > 1 else np.array([[1.0]])
        h = bond(d0, d1, W, 1.0, 0.0, mu)
        xs = np.linspace(0, 1, 200001); e0 = np.min(mu*xs-(2/np.pi)*np.sin(np.pi*xs))
        for L in range(2, 7):
            if d**L > 1100: continue
            E = np.linalg.eigvalsh(chain_H(h, d, L, True))[0]
            val = L*(L*e0-E); bnd = 2*L*(1/np.sin(np.pi/L)-L/np.pi)
            worst = max(worst, val/bnd)
            assert val <= bnd+1e-7, (d0, d1, mu, L, val, bnd)
print(f"(B2): all cases pass; largest ratio L(L e0-E_L)/[2tL(csc(pi/L)-L/pi)] = {worst:.6f}; constant (4-8/pi) = {4-8/np.pi:.5f}")

# SU(3) swap control: H = sum P_{j,j+1}; M = 1; e0 = 1 - ln3 - pi/(3 sqrt 3) (Sutherland)
e0 = 1-np.log(3)-np.pi/(3*np.sqrt(3))
P = np.zeros((9, 9))
for a in range(3):
    for b in range(3): P[b*3+a, a*3+b] = 1
for L in (3, 6, 9):
    dim = 3**L; H = sp.csr_matrix((dim, dim))
    idx = np.arange(dim); digs = np.array([(idx//3**(L-1-k)) % 3 for k in range(L)])
    for j in range(L):
        a, b = j, (j+1) % L
        new = digs.copy(); new[[a, b]] = digs[[b, a]]
        tgt = sum(new[k]*3**(L-1-k) for k in range(L))
        H = H+sp.csr_matrix((np.ones(dim), (tgt, idx)), shape=(dim, dim))
    E = np.linalg.eigvalsh(H.toarray())[0] if dim < 2000 else spl.eigsh(H, k=1, which='SA')[0][0]
    print(f"SU(3) swap L={L}: E_L={E:.8f}  L(L e0-E_L)/M={L*(L*e0-E):.4f}  (limit pi A/6 = pi^2*2/9 = {2*np.pi**2/9:.4f}; (B2) constant 1.4535)")

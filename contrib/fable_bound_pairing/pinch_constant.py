"""Sharp constant of the site-1 pinching removal on C^d x C^d (numerical).

Maximizes r(h) = ||h - pinch_1(h)|| / ||h|| over Hermitian h on C^d (x) C^d
(site-2 factor included) by random restarts + Nelder-Mead. Expected sharp
value 2(1 - 1/d): attained already by h = (J - (d/2) 1) (x) 1 with J the
all-ones matrix. Only numpy and scipy."""
import numpy as np
from scipy.optimize import minimize


def pinch1(h, d):
    hh = h.reshape(d, d, d, d).copy()
    for a in range(d):
        for ap in range(d):
            if a != ap:
                hh[ap, :, a, :] = 0
    return hh.reshape(d * d, d * d)


def opnorm(X):
    return np.linalg.norm(X, 2)

rng = np.random.default_rng(7)


def herm(p, n):
    A = np.zeros((n, n), complex)
    iu = np.triu_indices(n, 1)
    k = len(iu[0])
    A[iu] = p[:k] + 1j * p[k:2 * k]
    A = A + A.conj().T
    A[np.diag_indices(n)] = p[2 * k:]
    return A


for d in (2, 3, 4):
    n = d * d
    k = n * (n - 1) // 2
    best = 0
    for trial in range(12 if d < 4 else 4):
        p0 = rng.normal(size=2 * k + n)
        f = lambda p: -opnorm((lambda h: h - pinch1(h, d))(herm(p, n))) / opnorm(herm(p, n))
        res = minimize(f, p0, method='Nelder-Mead', options={'maxiter': 20000, 'xatol': 1e-9, 'fatol': 1e-12})
        best = max(best, -res.fun)
    print(f"d = {d}: max ratio found {best:.5f}, candidate sharp value 2(1-1/d) = {2*(1-1/d):.5f}")

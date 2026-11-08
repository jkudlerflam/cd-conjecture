#!/usr/bin/env python3
"""C165: residual uniform on-site symmetries of a bond (J, B): minimise ||Q^T J Q - J||^2 + ||Q B - B||^2
over Q in SO(3) (unitary) and over Q in O(3) minus SO(3) (antiunitary), excluding Q near 1.
Returns the smallest residual found away from the identity (0 means a symmetry exists)."""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation as R

def resid(J, B, Q):
    return np.sum((Q.T @ J @ Q - J) ** 2) + np.sum((Q @ B - B) ** 2)

def check(J, B, nstart=60, seed=0):
    rng = np.random.default_rng(seed); out = {}
    for name, s in (('unitary', 1), ('antiunitary', -1)):
        best = np.inf
        for _ in range(nstart):
            x0 = rng.standard_normal(3)
            f = lambda x: resid(J, B, s * R.from_rotvec(x).as_matrix()) + (1.0 if (s == 1 and min(np.linalg.norm(x) % (2*np.pi), 2*np.pi - np.linalg.norm(x) % (2*np.pi)) < 0.3) else 0.0)
            r = minimize(f, x0, method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-16, maxiter=4000))
            best = min(best, r.fun)
        out[name] = best
    return out

if __name__ == '__main__':
    from classical import fam
    a = 2 * np.pi * 13 / 100
    for kw in (dict(g=0.3, hx=0.2), dict(g=0.3, hx=0.2, hz=0.1), dict(g=0.3, hx=0.2, hz=0.1, Dx=0.05), dict(g=0.3), dict()):
        J, B = fam(a=a, Dz=1.02, **kw)
        print(kw, check(J, B))

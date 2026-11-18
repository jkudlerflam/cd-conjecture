"""C200: true Bethe roots at general filling against Theorem 4.11 and Corollary 4.12.
sup_v L|N_x - N_0| (evaluated at the roots, both one-sided values) and eta L sup|sigma_L - s| on a grid.  Floating point."""
import numpy as np, json
from scipy.interpolate import CubicHermiteSpline
from bethe_ferro import Model, solve
from fftsolver import Kern, solve_Q2
from gen_trial import counting_function
rows = []
for eps in [1e-2, 1e-3, 1e-4, 1e-5]:
    eta = np.arccos(1-eps); m = Model(eta); K = Kern(eta)
    for L in [64, 256, 1024]:
        for n in [L//8, L//4, 3*L//8, L//2-L//16, L//2-2]:
            nu = n/L
            R = solve_Q2(K, nu, hfac=10.0, withZ=False); Q = R.Q
            u, Z, s = counting_function(K, m, R, 25.0); sp = CubicHermiteSpline(u, Z, s)
            G, x, info = solve(L, n, eta)
            tN = np.clip(L*(sp(x)+nu/2), 0, n); j = np.arange(1, n+1)
            supc = max(np.max(np.abs(j-tN)), np.max(np.abs(j-1-tN)))
            X = np.linspace(-(Q+3), Q+3, 16001)
            T = lambda w: np.sin(2*eta)/(np.pi*(np.cosh(2*w)-np.cos(2*eta)))
            sig = np.sin(eta)/(np.pi*(np.cosh(2*X)+np.cos(eta)))
            for a in range(0, n, 256):
                sig += T(X[:, None]-x[None, a:a+256]).sum(axis=1)/L
            dens = eta*L*np.max(np.abs(sig-sp(X, 1)))
            rows.append(dict(eps=eps, L=L, n=n, lamQ=eta*L*R.r[-1], supcount=float(supc), dens=float(dens)))
            print(f'1+D={eps:g} L={L} n={n} lamQ={eta*L*R.r[-1]:.4f}: sup L|N_x-N_0| = {supc:.3f}   eta L sup|sigma-s| = {dens:.3f}', flush=True)
json.dump(rows, open('check_gen_roots.json', 'w'), indent=1)
print('max supcount', max(r['supcount'] for r in rows), 'max dens', max(r['dens'] for r in rows))

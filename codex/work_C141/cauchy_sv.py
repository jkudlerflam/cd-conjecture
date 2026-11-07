"""Singular-value decay of second-order port matrices (PROOF.md section 6.2).
T_ij = x_i y_j /(E_i + F_j), E_i, F_j in [a, b] (intermediate energies of two adjacent cells).
Cited bound (Beckermann-Townsend 2017, Cauchy matrices, rank-one displacement):
  s_{k+1}/s_1 <= Z_k <= 4 exp(-pi^2 k / log(4 b / a)).
We test it on box-like spectra E_i = a + (b - a) ((i-1)/(n-1))^2 and on log-uniform spectra,
with random and with end-visibility weights x_i ~ i/n^{3/2}.
"""
import numpy as np
rng = np.random.default_rng(1)
worst = 0.0
for a, b, n in [(1e-2, 1, 60), (1e-3, 1, 120), (1e-4, 1, 200), (1e-6, 1, 300)]:
    for kind in ['box', 'log']:
        if kind == 'box': E = a + (b-a)*(np.arange(n)/(n-1))**2
        else: E = np.exp(np.linspace(np.log(a), np.log(b), n))
        F = E[rng.permutation(n)]
        for wts in ['random', 'visibility']:
            x = rng.normal(size=n) if wts == 'random' else (np.arange(1, n+1)/n**1.5)
            y = rng.normal(size=n) if wts == 'random' else (np.arange(1, n+1)/n**1.5)
            T = np.outer(x, y)/(E[:, None]+F[None, :])
            s = np.linalg.svd(T, compute_uv=False)
            ks = np.arange(len(s))
            bound = 4*np.exp(-np.pi**2*ks/np.log(4*b/a))
            ratio = s/s[0]
            mask = ratio > 1e-13          # stay above double-precision floor
            viol = np.max(ratio[mask]/bound[mask])
            worst = max(worst, viol)
            keff = int(np.sum(ratio > 1e-8))
            print(f"a/b={a/b:.0e} n={n} {kind:3s} {wts:10s}: #(s/s1>1e-8)={keff:3d}  "
                  f"pred<= {1+np.log(4e8)*np.log(4*b/a)/np.pi**2:5.1f}  max ratio/bound={viol:.3f}")
print("worst ratio/bound:", worst)

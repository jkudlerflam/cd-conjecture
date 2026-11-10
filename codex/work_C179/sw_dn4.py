"""From the exact effective Hamiltonian of sw_exact.py (L = 11 open chain, w = i t2, Delta = 1):
(i) the largest Delta N = 4 matrix element, predicted O(Delta^6 / V^5) (no connected lower-order path into P);
(ii) the largest Delta N = 2 element between states whose particle-position sums differ by an odd number
     (U-odd under inversion x staggered rotation), predicted zero to all orders on Jx + Jy = 0;
(iii) the largest Delta N = 2 element overall, predicted O(Delta^3/V^2)."""
import numpy as np
from sw_exact import heff
L = 11
for V in (20.0, 40.0, 80.0):
    t2 = 1 / (2 * V); w = 1j * t2; mu = -2.2 * t2 + 3 / V
    He, pos, gap, smin = heff(L, 1.0, V, w, mu)
    S = list(pos.keys())
    occ = lambda s: [j for j in range(L) if (s >> (L - 1 - j)) & 1]
    N = np.array([len(occ(s)) for s in S]); X = np.array([sum(occ(s)) for s in S])
    dN = N[:, None] - N[None, :]
    A = np.abs(He)
    m4 = A[np.abs(dN) == 4].max()
    odd = (np.abs(dN) == 2) & ((X[:, None] - X[None, :]) % 2 == 1)
    even = (np.abs(dN) == 2) & ((X[:, None] - X[None, :]) % 2 == 0)
    print(f"V={V:5.0f}: max|dN=4| = {m4:.3e} (Delta^6/V^5 = {V**-5:.3e}); max|dN=2, odd position-sum change| = {A[odd].max():.3e} (Delta^3/V^2 = {V**-2:.3e}); "
          f"max|dN=2, even position-sum change| = {A[even].max():.3e}")

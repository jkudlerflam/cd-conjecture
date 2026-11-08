"""C162. Case B probe in C161 notation W = -h n_x + (kappa/2) n_z^2 + (mu/2) n_y^2 with mu < 0 (easy in-plane y):
flat in-plane background (m_u^2 = h + mu small), stiff out-of-plane background (m_w^2 = h + kappa).
For |mu| -> h the planar onset is at h = 1 exactly (sqrt(2 dV) = sqrt(h)(1 - cos phi)), and the core at phi = pi has
U = 4 sqrt(h) - 5 h + kappa, so the core is tilt-soft for kappa < 1.  Optional off-axis field hy breaks x -> -x, n_y -> -n_y."""
import sys
import numpy as np
from kfun import KFun, planar, tilt_lowest

hy = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0
print('mu     kappa  hy     s_c     m_u^2   m_w^2   lam0     glob')
for mu in (-0.7, -0.8, -0.9):
    for kappa in (0.0, 0.2, 0.4, 0.6):
        F = KFun(-kappa, mu - kappa, 1.0, hy, 1.0)
        pb, glob = F.background()
        mu2, mw2 = F.masses(pb)
        Egp, xs, p = planar(F, pb, X=min(40.0, 10 / np.sqrt(mu2)), dx=0.1)
        sc = (2 * np.pi / Egp) ** 2
        Fs = KFun(-kappa, mu - kappa, 1.0, hy, sc)
        Egp2, xs2, p2 = planar(Fs, pb, X=min(40.0, 10 / np.sqrt(mu2 * sc)), dx=0.1)
        lam0, _ = tilt_lowest(Fs, pb, xs2, p2)
        print(f'{mu:+.2f}  {kappa:+.2f}  {hy:+.2f}  {sc:.4f}  {mu2*sc:.4f}  {mw2*sc:.4f}  {lam0:+.4f}  {glob}')
        sys.stdout.flush()

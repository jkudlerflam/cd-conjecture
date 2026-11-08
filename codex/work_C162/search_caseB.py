"""C162. Search the K family (kfun.py) for shapes with
   (i) in-plane background that is the global minimum of W on S^2 and Z2-symmetric (n_z = 0),
   (ii) slow in-plane tail: m_w^2 > m_u^2 (Case B, out-of-plane tail decays faster),
   (iii) planar soliton tilt-unstable at its onset: lowest out-of-plane Hessian eigenvalue lam0 < 0.
The planar soliton's DM energy is -2 pi (topological), so its gradient+potential part scales as sqrt(s) and the
planar onset scale is s_c = (2 pi / E_gp(s=1))^2 exactly."""
import sys
import numpy as np
from kfun import KFun, planar, tilt_lowest

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
nsamp = int(sys.argv[2]) if len(sys.argv) > 2 else 150
print('a       b       hx      hy      s_c      mu2(s_c)  mw2(s_c)  ratio   lam0(s_c)')
for it in range(nsamp):
    a, b = rng.uniform(-3, 3, 2); hx = rng.uniform(0, 2); hy = rng.uniform(-2, 2)
    F = KFun(a, b, hx, hy, 1.0)
    pb, glob = F.background()
    if not glob:
        continue
    mu2, mw2 = F.masses(pb)
    if mu2 <= 0.02 or mw2 <= mu2:
        continue
    Egp, xs, p = planar(F, pb, X=min(40.0, 10.0 / np.sqrt(mu2)), dx=0.1)
    sc = (2 * np.pi / Egp) ** 2
    Fs = KFun(a, b, hx, hy, sc)
    # rescaled profile: phi_s(x) = phi_1(sqrt(s) x); recompute directly at s_c
    Egp2, xs2, p2 = planar(Fs, pb, X=min(40.0, 10.0 / np.sqrt(mu2 * sc)), dx=0.1)
    lam0, _ = tilt_lowest(Fs, pb, xs2, p2)
    tag = '  <-- tilt-unstable' if lam0 < 0 else ''
    print(f'{a:+.3f} {b:+.3f} {hx:+.3f} {hy:+.3f}  {sc:7.4f}  {mu2*sc:8.4f}  {mw2*sc:8.4f}  {mw2/mu2:6.2f}  {lam0:+.4f}  (E_gp-2pi={Egp2-2*np.pi:+.1e}){tag}')
    sys.stdout.flush()

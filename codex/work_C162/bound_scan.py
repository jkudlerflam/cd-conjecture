"""C162. Random K shapes (a, b, hx, hy) on the unit sphere, background = Z2-symmetric global minimum.
At the planar onset scale s_c record lam0 (tilt eigenvalue of the planar soliton), m_u^2, m_w^2, and
P_min = min_x [2 phi' + 2 (L(phi) - L_bg)], the term that enters the exact identity
    H_perp = H_par + (m_w^2 - m_u^2) + 2 phi' + 2 (L(phi) - L_bg)        (PROOF section 4.2),
so lam0 >= m_w^2 - m_u^2 whenever P_min >= 0.  Reports how often tilt instability occurs with m_w >= m_u."""
import sys, time
import numpy as np
from kfun import KFun, planar, tilt_lowest

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
nsamp = int(sys.argv[2]) if len(sys.argv) > 2 else 300
t0 = time.time(); rows = []
for it in range(nsamp):
    q = rng.normal(size=4); q /= np.linalg.norm(q); q[2] = abs(q[2])
    a, b, hx, hy = q
    F = KFun(a, b, hx, hy, 1.0)
    pb, glob = F.background()
    mu2, mw2 = F.masses(pb)
    if not glob or mu2 < 0.02 or mw2 < 0.005:
        continue
    Egp, xs, p = planar(F, pb, X=min(30.0, 8 / np.sqrt(mu2)), dx=0.12)
    sc = (2 * np.pi / Egp) ** 2
    Fs = KFun(a, b, hx, hy, sc)
    Egp2, xs2, p2 = planar(Fs, pb, X=min(30.0, 8 / np.sqrt(mu2 * sc)), dx=0.1)
    lam0, _ = tilt_lowest(Fs, pb, xs2, p2)
    dx = xs2[1] - xs2[0]
    d = np.diff(np.r_[pb, p2, pb + 2 * np.pi]) / dx; dp = 0.5 * (d[:-1] + d[1:])
    L = sc * (hx * np.cos(p2) + hy * np.sin(p2)); Lb = sc * (hx * np.cos(pb) + hy * np.sin(pb))
    Pmin = np.min(2 * dp + 2 * (L - Lb))
    rows.append((mu2 * sc, mw2 * sc, lam0, Pmin, a, b, hx, hy))
rows = np.array(rows)
mu2, mw2, lam0, Pmin = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]
print(f'{len(rows)} admissible shapes, {time.time()-t0:.0f} s')
print(f'tilt-unstable (lam0 < 0): {np.sum(lam0 < 0)};  of these with m_w^2 >= m_u^2: {np.sum((lam0 < 0) & (mw2 >= mu2))}')
print(f'min over shapes of lam0 - (m_w^2 - m_u^2): {np.min(lam0 - (mw2 - mu2)):+.4f}')
print(f'shapes with P_min < 0: {np.sum(Pmin < -1e-6)};  among them lam0 - (m_w^2 - m_u^2) min = '
      f'{np.min((lam0 - (mw2 - mu2))[Pmin < -1e-6]) if np.any(Pmin < -1e-6) else float("nan"):+.4f}')
print('tilt-unstable shapes (m_u^2, m_w^2, lam0, P_min, a, b, hx, hy):')
for r in rows[lam0 < 0][:12]:
    print('  ' + '  '.join(f'{v:+.4f}' for v in r))

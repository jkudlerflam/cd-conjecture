"""C162. Targeted search: minimize lam0(s_c) (planar-soliton tilt eigenvalue at the planar onset) over the K shapes
(a, b, hx, hy) subject to a Z2-symmetric global-minimum background and m_w^2 >= r * m_u^2 (r = argv[1]).
Penalty method, Nelder-Mead from several starts."""
import sys
import numpy as np
from scipy.optimize import minimize
from kfun import KFun, planar, tilt_lowest

r = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
def obj(q, verbose=False):
    q = np.asarray(q) / np.linalg.norm(q)
    a, b, hx, hy = q
    F = KFun(a, b, hx, hy, 1.0)
    pb, glob = F.background()
    mu2, mw2 = F.masses(pb)
    pen = 0.0
    if not glob: pen += 5.0
    if mu2 < 0.03: pen += 5.0 * (0.03 - mu2) / 0.03 + 1
    if mw2 < r * mu2: pen += 10.0 + 10.0 * (r * mu2 - mw2) / max(mu2, 1e-3)
    if pen > 0: return pen + 10
    try:
        Egp, xs, p = planar(F, pb, X=min(30.0, 8 / np.sqrt(max(mu2, 1e-3))), dx=0.12)
        sc = (2 * np.pi / Egp) ** 2
        Fs = KFun(a, b, hx, hy, sc)
        Egp2, xs2, p2 = planar(Fs, pb, X=min(30.0, 8 / np.sqrt(max(mu2 * sc, 1e-3))), dx=0.12)
        lam0, _ = tilt_lowest(Fs, pb, xs2, p2)
    except Exception:
        return 20.0
    val = lam0      # at the onset scale; < 0 means tilt-unstable
    if verbose:
        print(f'  a={a:+.3f} b={b:+.3f} hx={hx:+.3f} hy={hy:+.3f}  s_c={sc:.4f} m_u^2={mu2*sc:.4f} m_w^2={mw2*sc:.4f} '
              f'ratio={mw2/mu2:.2f} lam0={lam0:+.4f} glob={glob}')
    return val

starts = [(-0.3, 2.0, 0.5, 0.5), (-0.5, 2.0, 1.0, 1.0), (0.5, 3.0, 2.0, 1.5), (-1.0, 1.0, 1.0, -0.5), (0.0, -0.9, 1.0, 0.2)]
for q0 in starts:
    res = minimize(obj, q0, method='Nelder-Mead', options=dict(maxfev=160, xatol=1e-3, fatol=1e-4))
    print(f'start {q0}: best lam0 (+ penalty if infeasible) = {res.fun:+.4f}')
    obj(res.x, verbose=True)
    sys.stdout.flush()

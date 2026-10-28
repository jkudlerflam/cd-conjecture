"""Independent check of the band analysis at small n.  Recomputes T, mu at the
schedule point read from results/n{n}.json, rescales to O(1), and finds Fermi
points and centres by brute force on a dense k grid with ycell.bands (double
precision, no whitening, no Newton).  Compares with the high-precision values.
Usage: python3 -B band_check.py [n ...]  (default 2 3 4 5).  Writes band_check.json."""
import sys, os, json
import numpy as np
from decimal import Decimal
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'routeA_spatial'))
import ycell                                   # noqa: E402
import hpmat                                   # noqa: E402
from leading import leading, default_params, beta   # noqa: E402
from run_schedule import prec_for              # noqa: E402

out = []
for n in [int(a) for a in sys.argv[1:]] or [2, 3, 4, 5]:
    R = json.load(open(os.path.join(HERE, 'results', 'n%d.json' % n)))
    run = R['runs']['C33_kappa1e-3']
    hpmat.setprec(prec_for(n))
    P = default_params(n)
    b = beta(n, P['a'], P['Dx'])
    P['td'] = b * Decimal(2) ** (-run['p'])
    eps = 2 * b * b * Decimal(2) ** (-run['f'])
    L = leading(n, P)
    sc = Decimal(L['T'].absmax())
    T = L['T'].rmul(1 / sc).tocomplex()
    mu = L['mu'].rmul(1 / sc).tocomplex() + float(eps / sc) * np.eye(n)
    bs = ycell.bands(T, mu, nk=2000001)
    # high-precision values
    hp_c = [float(x['centre_k']) for x in run['bands']]
    hp_kF = [float(x['kF']) for x in run['bands']]
    # ycell orders bands bottom-up, as here.  Its own centre is a grid average
    # (resolution 2 pi / nk); use the midpoint of the occupied arc between the
    # two interpolated roots instead, and ycell's centre only to pick the arc.
    grid_c, grid_kF = [], []
    for x in bs:
        a, c = sorted(x['roots'])
        if x['centre'] > a and x['centre'] < c:
            mid, arc = (a + c) / 2, c - a
        else:
            mid, arc = (a + c) / 2 + np.pi, 2 * np.pi - (c - a)
            mid = (mid + np.pi) % (2 * np.pi) - np.pi
        grid_c.append(mid)
        grid_kF.append(arc / 2)
    rec = dict(n=n, roots_per_band=[len(x['roots']) for x in bs],
               max_centre_diff=float(max(abs(a - c) for a, c in zip(hp_c, grid_c))),
               max_kF_diff=float(max(abs(a - c) for a, c in zip(hp_kF, grid_kF))),
               min_centre_sep_grid=float(min(abs(grid_c[i] - grid_c[j]) for i in range(n) for j in range(i + 1, n))) if n > 1 else None,
               min_centre_sep_hp=float(min(abs(hp_c[i] - hp_c[j]) for i in range(n) for j in range(i + 1, n))) if n > 1 else None,
               narrow_over_wide=float(Decimal(run['lambda_min_Rs']) / Decimal(run['lambda_max_Rs'])))
    print(json.dumps(rec), flush=True)
    out.append(rec)
json.dump(out, open(os.path.join(HERE, 'band_check.json'), 'w'), indent=1)

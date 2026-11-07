"""Smallest majorant constants at m = infinity.
For a given (S, V') the smallest admissible K is sup over (L, beta) of log Z~(L) - S F(L/(V' beta)).
Sources: (a) small rings L <= 14, beta <= 8, from data/infinite_m.json (exact ring spectra plus extrapolated eps0);
(b) L = 2 for 8 <= beta <= 1e4 from the 2-ring spectrum (exactly 2h) and the low-temperature form
    eps0(beta) = beta e0 - pi c/(6 v beta) (+ O(beta^-2)), where only c v = A enters; the error of this form at
    beta = 8 is measured against (a) and reported;  (c) XX: exact for all (L, beta) (xx_exact.npz).
Prints K_min(V') for S = 2 (and S = 1) on a grid of V'/M and V'_min/M for K = log 2, 1, 2."""
import os, json
import numpy as np
from scipy.special import logsumexp
from majorant import F
from models import get, MODELS

HERE = os.path.dirname(os.path.abspath(__file__))


def L2_large_beta(name, betas):
    h, mode, M, A, label = get(name)
    z = json.load(open(os.path.join(HERE, 'data', 'zeroT.json')))[name]
    e0 = z['e0']
    lam = np.linalg.eigvalsh(2 * h)        # the L = 2 ring is H_2 = h_12 + h_21 = 2h (h symmetric)
    out = []
    for b in betas:
        eps0 = b * e0 - np.pi * A / (6 * b)
        out.append(logsumexp(-b * lam) + 2 * eps0)
    return np.array(out)


def data_points(name):
    """list of (L, beta, log Z~) at m = infinity"""
    pts = []
    im = json.load(open(os.path.join(HERE, 'data', 'infinite_m.json')))[name]
    for b, r in im.items():
        for L, g in r['logZt'].items():
            pts.append((int(L), float(b), g))
    bl = np.geomspace(8, 1e4, 400)
    for b, g in zip(bl, L2_large_beta(name, bl)):
        pts.append((2, float(b), float(g)))
    return np.array(pts)


def check_L2_asymptotic(name):
    im = json.load(open(os.path.join(HERE, 'data', 'infinite_m.json')))[name]
    exact8 = im['8']['logZt']['2']
    return exact8 - L2_large_beta(name, [8.0])[0]


def Kmin(pts, S, V):
    L, b, g = pts.T
    return np.max(g - S * F(L / (V * b)))


if __name__ == '__main__':
    res = {}
    for name in MODELS:
        h, mode, M, A, label = get(name)
        pts = data_points(name)
        z = json.load(open(os.path.join(HERE, 'data', 'zeroT.json')))[name]
        row = dict(asymptotic_check_L2_beta8=check_L2_asymptotic(name))
        for S in (1, 2):
            V0 = 12 * z['Dmax'] / (np.pi * S)
            Vs = M * np.linspace(1.0, 20, 2000)
            Vs = Vs[Vs > V0 * 1.0005]
            Ks = np.array([Kmin(pts, S, V) for V in Vs])
            for K in (np.log(2), 1.0, 2.0):
                ok = Vs[Ks <= K]
                row['S%d_Vmin_over_M_K%.2f' % (S, K)] = float(ok.min() / M) if ok.size else None
            row['S%d_V0_over_M' % S] = V0 / M
        res[name] = row
        print(name, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()})
    json.dump(res, open(os.path.join(HERE, 'data', 'minimal_constants.json'), 'w'), indent=1)

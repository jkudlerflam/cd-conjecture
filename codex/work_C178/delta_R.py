"""Effective pseudospin anisotropy from the ratio R = D_neel / D_ps of the twist-relaxed spectrum.

For an N-site XXZ ring with twist Phi, R(Delta, Phi) = E(Sz = 0, q = pi) / E(Sz = 1) - both relative to the ground
state - is 1 at the SU(2) point for Phi = 0, above 1 for easy-plane, below 1 for easy-axis (xxz_fit.features). The
residual pseudospin twist is taken as Phi_s = (k_b - k_a) L mod 2 pi with the bare valley minima (PROOF 2). Delta_R is
the anisotropy that reproduces the measured R at that twist (linear interpolation on a grid). Also printed: the
log-slope of D_ps between consecutive even-N sizes at the same (|w|, n), the pre-registered P1 statistic.
Usage: python3 delta_R.py files.jsonl"""
import sys, json, numpy as np
from collections import defaultdict
from fractions import Fraction
import analyze2 as A, xxz_fit as X

grid = np.linspace(-0.9, 3.0, 79)
cache = {}


def Rcurve(N, Phi):
    key = (N, round(Phi, 6))
    if key not in cache:
        cache[key] = np.array([X.features(N, d, Phi)[f'k{N // 2}'] / X.features(N, d, Phi)['D_ps'] for d in grid])
    return cache[key]


rows = defaultdict(list)
for fn in sys.argv[1:]:
    for line in open(fn):
        r = json.loads(line)
        if r['model'] != 'cusp' or r['N'] % 2 or r['N'] > 10:
            continue
        o = A.summarize(r, verbose=False)
        if o.get('R') is None:
            continue
        ks, es, cs = A.valley_cos(r['phi'], r['wabs'])
        Phi = (abs(ks[0] - ks[1]) * r['L']) % (2 * np.pi)
        Rc = Rcurve(r['N'], Phi)
        # R decreases monotonically with Delta on the grid; invert
        o_R = o['R']
        if o_R > Rc[0]: dR = f"< {grid[0]:.1f}"
        elif o_R < Rc[-1]: dR = f"> {grid[-1]:.1f}"
        else:
            i = np.nonzero(Rc >= o_R)[0][-1]
            j = min(i + 1, len(grid) - 1)
            t = 0 if Rc[i] == Rc[j] else (Rc[i] - o_R) / (Rc[i] - Rc[j])
            dR = f"{grid[i] + t * (grid[j] - grid[i]):.2f}"
        R1 = Rcurve(r['N'], Phi)[np.argmin(abs(grid - 1.0))]
        rows[(round(r['phi'], 4), round(r['wabs'], 4), Fraction(r['N'], r['L']))].append(
            (r['L'], r['N'], o['D_ps'], o['D_neel'], o['R'], Phi / np.pi, R1, dR))
for k in sorted(rows, key=lambda k: (k[0], k[1], float(k[2]))):
    print(f"\nphi={k[0]} |w|={k[1]} n={k[2]}")
    prev = None
    uniq = {}
    for row in rows[k]: uniq[(row[0], row[1])] = row
    for (L, N, Dps, Dn, R, Phi, R1, dR) in sorted(uniq.values()):
        sl = f"{np.log(Dps / prev[1]) / np.log(L / prev[0]):+.2f}" if prev else "  -  "
        print(f"   L={L:3d} N={N:2d}: D_ps={Dps:.4f} L*D_ps={L*Dps:.3f} slope={sl}  D_neel={Dn:.4f} R={R:.3f}  Phi_s={Phi:.3f}pi  R(SU2 at this Phi)={R1:.3f}  Delta_R={dR}")
        prev = (L, Dps)

#!/usr/bin/env python3
"""C177: tune D (so that the exact one-particle band has |w|/|t2| = 1) and h_z (so that the dressed one-particle
chemical potential mu_eff = -c0/|t2| equals a target) for the full chain at V/Delta = r, phi = pi/2 stratum.
c0(hz) is measured at two values of hz far below the band (one-particle sector isolated) and extrapolated
linearly to the target (the dependence is linear up to O(t2 Delta^2/V^2) dressing changes, checked by a third point).
Ring L = 12 (L = 10, 12, 14 agree to 1e-6 relative at r = 8, band.out)."""
import sys, json
import numpy as np
import band as B

def measure(r, Dscale, mu_eff, L=12, Jxpy=0.0):
    cp = B.couplings(r, mu_eff=mu_eff, Dscale=Dscale, Jxpy=Jxpy)
    Ev, eps = B.band(B.bond(**cp), L)
    f, res = B.fit(eps, L)
    t2 = -f['b2'] / 2
    return cp, f, t2, abs(f['a1'] / 2) / t2, eps

def tune(r, target_mu=-2.2, L=12):
    Ds = 1.0
    for it in range(3):
        cp, f, t2, ratio, _ = measure(r, Ds, -12.0, L)
        Ds = Ds / ratio
    # c0 is -mu_dressed; get hz giving mu_dressed/|t2| = target
    pts = []
    for mu_lo in (-12.0, -9.0, -6.0):
        cp, f, t2, ratio, eps = measure(r, Ds, mu_lo, L)
        pts.append((cp['hz'], f['c0'], t2, ratio))
    hz = np.array([p[0] for p in pts]); c0 = np.array([p[1] for p in pts])
    a, b = np.polyfit(hz, c0, 1)
    quad = np.polyfit(hz, c0, 2)[0]
    t2 = pts[0][2]
    hz_t = (-target_mu * t2 - b) / a
    cp = B.couplings(r, mu_eff=target_mu, Dscale=Ds)
    cp['hz'] = float(hz_t)
    return dict(r=r, Dscale=Ds, t2=t2, t2_over_lo=t2 * 2 * r, ratio=pts[0][3], dc0_dhz=a, quad=quad,
                lin_resid=float(np.max(np.abs(np.polyval([a, b], hz) - c0))), couplings=cp, target_mu=target_mu,
                hz_leading=B.couplings(r, mu_eff=target_mu)['hz'])

if __name__ == '__main__':
    out = []
    for r in [float(x) for x in sys.argv[1].split(',')]:
        d = tune(r); out.append(d); print(json.dumps(d), flush=True)
    if len(sys.argv) > 2:
        json.dump(out, open(sys.argv[2], 'w'), indent=1)

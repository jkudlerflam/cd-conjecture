"""Bias of the local demodulation filter of ana_quartet.py on the exponents of REPORT 9.217 (2026-11-25, Claude Opus 5.5).

Independent of contrib/quartet_typhon_20261125: it reuses only k0() and quad() of ana_quartet.py and the committed data.
For each correlator it (1) reproduces the exponent quoted in tables/quartet_k0.txt, a straight-line fit of the filtered
correlator over r = 32 to 56, (2) passes a pure power r^-p through the same filter on the same distance grid (same
starting distance, same wavevector, same centres), and (3) finds the p whose filtered version has the quoted slope.
Step (3) assumes a pure power law over the whole support of the filter (from the operator's own starting distance on),
which the data need not obey, so the result is a corrected reading and not a measured dimension.  The spread over
filter widths Wd is printed because a pure power law would give the same corrected p at every width.

usage: python check_filter_bias.py            (from anywhere; needs numpy and scipy)
"""
import gzip, json, os, sys
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__)); sys.argv = sys.argv[:1]; sys.path.insert(0, HERE)
import ana_quartet as A

D = os.path.join(HERE, "..", "meas", "quartet")
K0 = A.k0
WIDTHS = (32, 40, 48, 56)


def use_width(Wd):
    A.k0 = lambda r, v, q0=0.0: K0(r, v, q0, Wd=Wd)


def slope(r, g, a=32, b=56):
    s = (r >= a) & (r <= b)
    same = np.all(np.sign(g[s]) == np.sign(g[s][0]))
    return (-np.polyfit(np.log(r[s]), np.log(np.abs(g[s])), 1)[0]) if same else np.nan


def filtered_power_slope(p, r0, n, Wd):
    r = np.arange(n) + r0
    rc, c = K0(r, r.astype(float) ** -p, 0.0, Wd=Wd)
    return slope(rc, c.real)


def corrected(y, r0, n, Wd):
    if not np.isfinite(y):
        return np.nan
    f = lambda p: filtered_power_slope(p, r0, n, Wd) - y
    return brentq(f, 0.2, 7.0) if f(0.2) < 0 < f(7.0) else np.nan


def row(m, blk, part, Wd):
    use_width(Wd)
    i0 = str(m["refs"][0]); r0 = blk[i0]["j0"] - m["refs"][0]
    n = min(len(blk[str(i)]["md"]) for i in m["refs"])
    r, big, small, _ = A.quad(blk, m["refs"])
    y = slope(r, big if part == "pair" else small)
    return r0, y, corrected(y, r0, n, Wd)


print("calibration (start r0 = 10, Wd = 48):", "  ".join("p=%.2f -> %.3f" % (p, (setattr(A, "KAPPA", 0.2077 * np.pi), filtered_power_slope(p, 10, 190, 48))[1]) for p in (1.0, 3.0, 4.0)))
for fn in ("T1_L256_chi384", "T1_L256", "T1_L384_chi512", "T1_L384"):
    m = json.load(gzip.open(os.path.join(D, fn + ".json.gz")))
    A.KAPPA = A.fit_kappa(m["P"]["2"], m["refs"][0])
    print("== %s  L=%d chi=%d" % (fn, m["L"], m["chi"]))
    for part, key, blk in (("pair", "2", m["P"]["2"]), ("pair", "3", m["P"]["3"]), ("quartet", "2_5", m["Q"]["2_5"]), ("quartet", "3_6", m["Q"]["3_6"])):
        out = []
        for Wd in WIDTHS:
            r0, y, p = row(m, blk, part, Wd)
            out.append("%d: %s -> %s" % (Wd, "%.2f" % y if np.isfinite(y) else "sign change", "%.2f" % p if np.isfinite(p) else "--"))
        print("   %-7s %-3s r0=%-2d  Wd: quoted slope -> corrected p    %s" % (part, key, r0, "    ".join(out)))

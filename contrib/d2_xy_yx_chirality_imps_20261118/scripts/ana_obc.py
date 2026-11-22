"""Open-chain c estimators from out_obc/*.npz.
usage: python ana_obc.py [substr ...] [Lmin=96] [win=0.125]

Per run (lab without _L): for each L the chi ladder (S_mid, trunc, conv), then at the largest chi
  c_prof : linear fit of S(l) against X(l) = (1/6) ln[(2L/pi) sin(pi l/L)] over L/8 <= l <= 7L/8
  Sbar   : mean of S(l) over |l - L/2| <= win*L,  Xbar the same mean of X
and across L:  c_L = slope of Sbar against Xbar for L >= Lmin (with error), and pairwise local slopes.
The chi-convergence of Sbar is shown as dSbar = Sbar(chi_max) - Sbar(previous chi)."""
import sys, glob, re
import numpy as np

args = [a for a in sys.argv[1:] if "=" not in a]
kw = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
Lmin = int(kw.get("Lmin", 96)); win = float(kw.get("win", 0.125))
runs = {}
for fn in glob.glob("out_obc/*.npz"):
    if args and not all(a in fn for a in args):
        continue
    m = re.match(r"out_obc/(.*)_L(\d+)_p(\d)_chi(\d+)\.npz", fn)
    if not m:
        continue
    runs.setdefault((m.group(1), int(m.group(3))), {}).setdefault(int(m.group(2)), {})[int(m.group(4))] = np.load(fn)


def Xl(L):
    l = np.arange(1, L)
    return l, np.log(2 * L / np.pi * np.sin(np.pi * l / L)) / 6


for (lab, par), byL in sorted(runs.items()):
    print("==", lab, "par=%d" % par)
    rows = []
    for L in sorted(byL):
        lad = byL[L]; chis = sorted(lad)
        l, X = Xl(L)
        sel = (l >= L / 8) & (l <= 7 * L / 8); mid = np.abs(l - L / 2) <= win * L
        sb = []
        for c in chis:
            S = lad[c]["S"]
            sb.append(S[mid].mean())
        d = lad[chis[-1]]; S = d["S"]
        cp = np.polyfit(X[sel], S[sel], 1)[0]
        ladder = " ".join("%d:%.5f%s" % (c, lad[c]["S"][L // 2 - 1], "" if bool(lad[c]["conv"]) else "*") for c in chis)
        dsb = sb[-1] - sb[-2] if len(sb) > 1 else np.nan
        print("  L=%4d E/L=%.10f trunc=%.1e  Smid[chi]: %s  Sbar=%.5f dSbar(chi)=%+.1e  c_prof=%.3f" %
              (L, float(d["E"]) / L, float(d["trunc"]), ladder, sb[-1], dsb, cp))
        rows.append((L, X[mid].mean(), sb[-1]))
    rows = np.array(rows)
    if len(rows) >= 2:
        loc = ["%d-%d:%.3f" % (rows[i, 0], rows[i + 1, 0], (rows[i + 1, 2] - rows[i, 2]) / (rows[i + 1, 1] - rows[i, 1]))
               for i in range(len(rows) - 1)]
        print("  local c from Sbar:", " ".join(loc))
    s = rows[rows[:, 0] >= Lmin] if len(rows) else rows
    if len(s) >= 3:
        p, cov = np.polyfit(s[:, 1], s[:, 2], 1, cov=True)
        print("  c_L (L>=%d, n=%d) = %.3f +- %.3f" % (Lmin, len(s), p[0], np.sqrt(cov[0, 0])))

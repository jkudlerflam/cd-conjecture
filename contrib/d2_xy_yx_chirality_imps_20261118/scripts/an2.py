"""Tables from out3/*.npz (final VUMPS states).  xi excludes parity-cat partners (|lambda| > 1-1e-7).
usage: python an2.py [substr ...] [mode=c|m]   prints per-run chi ladder with local c_eff = 6 dS/dln xi and
local order-parameter exponent  dln|m| / dln xi."""
import sys, glob, re
import numpy as np

args = [a for a in sys.argv[1:] if "=" not in a]
kw = dict(a.split("=", 1) for a in sys.argv[1:] if "=" in a)
runs = {}
for fn in glob.glob("out3/*.npz"):
    if fn.endswith(".sq.npz") or (args and not all(a in fn for a in args)):
        continue
    m = re.match(r"out3/(.*)_chi(\d+)(_.*)\.npz", fn)
    if not m:
        continue
    lab, chi, tag = m.group(1), int(m.group(2)), m.group(3)
    if tag.endswith("pinned"):
        continue
    d = np.load(fn)
    tm0 = d["tm0"][1:] if len(d["tm0"]) > 1 else np.array([np.nan])
    tm1 = d["tm1"] if len(d["tm1"]) else np.array([np.nan])
    a0, a1 = np.abs(tm0), np.abs(tm1)
    good = np.concatenate([a0[a0 < 1 - 1e-7], a1[a1 < 1 - 1e-7]])
    xi = -2 / np.log(np.nanmax(good))
    xi0 = -2 / np.log(np.nanmax(a0[a0 < 1 - 1e-7])) if np.any(a0 < 1 - 1e-7) else np.nan
    xi1 = -2 / np.log(np.nanmax(a1[a1 < 1 - 1e-7])) if np.any(a1 < 1 - 1e-7) else np.nan
    cat = int(np.sum(a1 > 1 - 1e-7))
    k = float(d["kstag"])
    nxu = float(np.mean(d["nx"])) if "nx" in d.files else np.nan
    Z = float(np.mean(d["Z"]))
    import os
    sq = fn[:-4] + ".sq.npz"
    Sc = float(np.mean(np.load(sq)["Sconf"])) if os.path.exists(sq) else np.nan
    runs.setdefault((lab, tag), []).append((chi, float(d["E"]), float(np.mean(d["S"])), xi, xi0, xi1, k, nxu, Z, cat, Sc))
for key in sorted(runs):
    R = sorted(runs[key])
    print("==", key[0], key[1])
    prev = None
    for r in R:
        chi, E, S, xi, xi0, xi1, k, nxu, Z, cat, Sc = r
        ce = mexp = cc = c0 = ""
        if prev is not None:
            dl = np.log(xi / prev[3])
            if dl != 0:
                ce = "%.3f" % (6 * (S - prev[2]) / dl)
                cc = "%.3f" % (6 * (Sc - prev[10]) / dl)
            dl0 = np.log(xi0 / prev[4])
            if dl0 != 0:
                c0 = "%.3f" % (6 * (S - prev[2]) / dl0)
                if abs(k) > 0 and abs(prev[6]) > 0:
                    mexp = "%.3f" % (np.log(abs(k) / abs(prev[6])) / dl)
        print("  chi=%4d E=%.12f S=%.5f Sc=%.5f xi=%8.1f xi0=%8.1f xi1=%8.1f cat=%d m=%+.5f B=%+.5f Z=%.4f c_xi=%6s c_xi0=%6s c_conf=%6s dlnm/dlnxi=%6s"
              % (chi, E, S, Sc, xi, xi0, xi1, cat, k, nxu, Z, ce, c0, cc, mexp))
        prev = r
    sel = [r for r in R if r[0] >= int(kw.get("chimin", 64))]
    if len(sel) >= 3:
        x = np.log([r[4 if kw.get("xi0") else 3] for r in sel]); y = np.array([r[10 if kw.get("conf") else 2] for r in sel])
        if len(sel) >= 4:
            p, cov = np.polyfit(x, y, 1, cov=True); err = 6 * np.sqrt(cov[0, 0])
        else:
            p = np.polyfit(x, y, 1); err = np.nan
        print("  fit chi>=%s: c=%.3f +- %.3f   (n=%d)" % (kw.get("chimin", 64), 6 * p[0], err, len(sel)))

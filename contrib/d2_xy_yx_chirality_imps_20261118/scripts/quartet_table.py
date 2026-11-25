"""Table of the k = 0 quartet and pi-pair exponents from meas_quartet.py output (2026-11-24, Claude Opus 5.5).
usage: python quartet_table.py <dir with T1_L256.json ...>  (gunzip meas/quartet/*.json.gz first)     (uses ana_quartet.py for the demodulation and quadratures)"""
import sys, os, json, numpy as np
D = sys.argv[1]; sys.argv = sys.argv[:1]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import ana_quartet as A


def pw(r, g, a, b):
    s = (r >= a) & (r <= b)
    return -np.polyfit(np.log(r[s]), np.log(np.abs(g[s])), 1)[0]


def slopes(r, g, L):
    out = []
    for a in (36, 44, 52, 64, 80):
        b = a * 1.4
        if b > L / 3: break
        s = (r >= a) & (r <= b)
        ok = s.sum() >= 3 and np.all(np.sign(g[s]) == np.sign(g[s][0]))
        out.append("%5.2f" % -np.polyfit(np.log(r[s]), np.log(abs(g[s])), 1)[0] if ok else "   --")
    return " ".join(out)


for fn in ["T1_L256_chi384", "T1_L256", "T1_L384_chi512", "T1_L384", "T3_L256", "T3_L384"]:
    m = json.load(open(os.path.join(D, fn + ".json"))); L = m["L"]; refs = m["refs"]
    A.KAPPA = A.fit_kappa(m["P"]["2"], refs[0])
    print("== %s  L=%d chi=%d  kappa/pi=%.4f" % (fn, L, m["chi"], A.KAPPA / np.pi))
    for d in ("2", "3"):
        r, big, small, _ = A.quad(m["P"][d], refs)
        print("   pair d=%s      G(32)=%.1e  p[32,56]=%.2f  local(36,44,52,64,80)= %s" % (d, abs(big[np.argmin(abs(r - 32))]), pw(r, big, 32, 56), slopes(r, big, L)))
    for k in ("2_5", "3_6", "2_4", "2_1"):
        r, big, small, _ = A.quad(m["Q"][k], refs)
        print("   quartet %-4s  G(32)=%.1e  p[32,56]=%.2f  local(36,44,52,64,80)= %s    (neutral quadrature p[32,56]=%.2f)" % (
            k, abs(small[np.argmin(abs(r - 32))]), pw(r, small, 32, 56), slopes(r, small, L), pw(r, big, 32, 56)))

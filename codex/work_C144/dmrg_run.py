"""DMRG of the full qubit chain on open chains for h = XZ + alpha ZX + beta YY + extra terms.
Usage: python3 dmrg_run.py L chis(comma) out.jsonl  'label;alpha;beta;extra' ...   (extra: name=coef,... or '-')
Records energy, S(l), c fits (parity-resolved, central window), and bond expectations a_j = <X_j Z_j+1>,
b_j = <Z_j X_j+1>, y_j = <Y_j Y_j+1>, with their staggered parts at the centre."""
import sys, json, time
import numpy as np
from xzlib import term, SX, SY, SZ
from dmrg import dmrg, expect_bonds

def cfit(S, L, lo=0.25, hi=0.75):
    l = np.arange(1, L)
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    m = (l > lo * L) & (l < hi * L)
    out = []
    for par in (0, 1):
        mm = m & (l % 2 == par)
        out.append(float(6 * np.polyfit(x[mm], S[mm], 1)[0]))
    return out

L = int(sys.argv[1]); chis = tuple(int(c) for c in sys.argv[2].split(",")); fout = sys.argv[3]
for spec in sys.argv[4:]:
    label, al, be, extra = spec.split(";")
    h = term("xz") + float(al) * term("zx") + float(be) * term("yy")
    if extra != "-":
        for it in extra.split(","):
            n, c = it.split("="); h = h + float(c) * term(n)
    t0 = time.time()
    r = dmrg(h, L, chis=chis)
    ex = expect_bonds(r["M"], [(SX, SZ), (SZ, SX), (SY, SY)]).real
    a, b, y = ex
    mid = slice(L // 2 - 4, L // 2 + 4)
    stag = lambda v: float(np.mean(v[mid] * (-1.0) ** np.arange(L - 1)[mid]))
    rec = dict(label=label, L=L, alpha=float(al), beta=float(be), extra=extra, chis=list(chis), E=r["E"], E_per_site=r["E"] / L,
               trunc=r["hist"][-1][2], dE_last=r["hist"][-1][1] - r["hist"][-2][1] if len(r["hist"]) > 1 else None,
               c_fit=cfit(r["S"], L), c_fit_wide=cfit(r["S"], L, 0.15, 0.85), S_mid=float(r["S"][L // 2 - 1]),
               S=[round(float(s), 6) for s in r["S"]], a_mid=float(np.mean(a[mid])), b_mid=float(np.mean(b[mid])),
               a_stag=stag(a), b_stag=stag(b), y_stag=stag(y),
               a=[round(float(v), 6) for v in a], b=[round(float(v), 6) for v in b], sec=round(time.time() - t0, 1))
    with open(fout, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(f"{label} L={L} a={al} b={be} {extra}: E/L={rec['E_per_site']:.8f} trunc={rec['trunc']:.1e} dE={rec['dE_last']} "
          f"c_fit={np.round(rec['c_fit'],3).tolist()} wide={np.round(rec['c_fit_wide'],3).tolist()} S_mid={rec['S_mid']:.4f} "
          f"a_stag={rec['a_stag']:+.4f} b_stag={rec['b_stag']:+.4f} a_mid={rec['a_mid']:+.4f} b_mid={rec['b_mid']:+.4f} ({rec['sec']}s)", flush=True)

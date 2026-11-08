"""C154: open-chain DMRG on the cyclic stratum h = a XZ + b YX + g ZY (+ optional field) (Claude Opus agent, 2026-11-07).
Usage: python3 dmrg_cyc.py L chis(comma) out.jsonl a,b,g[,Bx,By,Bz] ...
Records E, entanglement S(l), bond energies e_j = <h_{j,j+1}>, the centre dimerization
D = mean over the central window of (-1)^j (e_j - e_{j+1})/2, and parity-resolved Calabrese-Cardy fits."""
import sys, json, time
import numpy as np
from xzlib import bond, SX, SY, SZ
from dmrg import dmrg, expect_bonds

def cfit(S, L, lo=0.25, hi=0.75):
    l = np.arange(1, L)
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    m = (l > lo * L) & (l < hi * L)
    return [float(6 * np.polyfit(x[m & (l % 2 == p)], S[m & (l % 2 == p)], 1)[0]) for p in (0, 1)]

L = int(sys.argv[1]); chis = tuple(int(c) for c in sys.argv[2].split(",")); fout = sys.argv[3]
for spec in sys.argv[4:]:
    v = [float(x) for x in spec.split(",")]
    a, b, g = v[:3]; B = v[3:6] if len(v) > 3 else [0, 0, 0]
    J = np.zeros((3, 3)); J[0, 2] = a; J[1, 0] = b; J[2, 1] = g
    h = bond(J, B)
    t0 = time.time()
    r = dmrg(h, L, chis=chis)
    xz, yx, zy = expect_bonds(r["M"], [(SX, SZ), (SY, SX), (SZ, SY)])
    e = (a * xz + b * yx + g * zy).real
    j = np.arange(L - 1); mid = (j > L // 4) & (j < 3 * L // 4 - 1)
    dim = float(np.mean(((-1.0) ** j[mid]) * (e[mid] - e[1:][mid[:-1]] if False else e[mid] - np.roll(e, -1)[mid]) / 2))
    Dc = float(abs(e[L // 2 - 1] - e[L // 2]) / 2)
    rec = dict(a=a, b=b, g=g, B=B, L=L, chis=list(chis), E=r["E"], E_per_site=r["E"] / L,
               trunc=r["hist"][-1][2], dE_last=(r["hist"][-1][1] - r["hist"][-2][1]) if len(r["hist"]) > 1 else None,
               c_fit=cfit(r["S"], L), c_fit_wide=cfit(r["S"], L, 0.15, 0.85), S=[round(float(s), 6) for s in r["S"]],
               e=[round(float(x), 7) for x in e], dimer_window=dim, dimer_centre=Dc, sec=round(time.time() - t0, 1))
    with open(fout, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(f"(a,b,g)=({a},{b},{g}) B={B} L={L}: E/L={rec['E_per_site']:.8f} trunc={rec['trunc']:.1e} dE={rec['dE_last']} "
          f"c_fit={np.round(rec['c_fit'],3).tolist()} wide={np.round(rec['c_fit_wide'],3).tolist()} S_mid={r['S'][L//2-1]:.4f} "
          f"dimer_centre={Dc:.5f} ({rec['sec']}s)", flush=True)

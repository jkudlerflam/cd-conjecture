"""C188: analysis of ring scans (ringk_run.py output). For each point and L: ground-state momentum m0, x = L(E-E0)
of the lowest levels, the SU(2)_1 identifications (triplet at m0+L/2+L/3, singlet at m0+L/2, current at m0+1),
v = x_c/2pi, and Casimir c from consecutive sizes. Usage: python3 analyze_scan.py file.jsonl [file2 ...]"""
import sys, json
import numpy as np
from collections import defaultdict
D = defaultdict(dict)
for fn in sys.argv[1:]:
    for l in open(fn):
        r = json.loads(l)
        D[(r['a'], r['b'], r['g'])].setdefault(r['L'], {})[r['m']] = r['levels']
def fold(m, L):
    m %= L
    return min(m, L - m)
def summary(sec, L):
    E0 = min(min(v) for v in sec.values())
    m0 = min(sec, key=lambda m: min(sec[m]))
    allv = sorted((x, m) for m, v in sec.items() for x in v)
    get = lambda m: sec.get(fold(m, L))
    tr = get(m0 + L // 2 + L // 3); pi = get(m0 + L // 2); cu = get(m0 + 1)
    if tr is None or pi is None or cu is None:
        return dict(E0=E0, m0=m0)
    t = tr[0]
    cand = [x for x in pi if abs(x - t) > 1e-7]
    s = cand[0] if cand else np.nan
    gap1 = L * (allv[1][0] - E0) if len(allv) > 1 else np.nan
    low = [(round(L * (x - E0), 3), m) for x, m in allv[:6]]
    return dict(E0=E0, m0=m0, x_t=L * (t - E0), x_s=L * (s - E0), x_c=L * (cu[0] - E0), gap1=gap1, low=low)
pts = sorted(D, key=lambda p: (p[1], -p[2]))
print("point (1,b,g) t=-g/b | L: m0 x_t x_s x_c gap1 | r_s r_c v | c_cas")
for p in pts:
    Ls = sorted(D[p])
    out = []
    S = {L: summary(D[p][L], L) for L in Ls}
    line = f"({p[0]},{p[1]},{p[2]}) t={-p[2]/p[1] if p[1] else 0:.2f} |"
    for L in Ls:
        s = S[L]
        if 'x_t' in s:
            line += f" L{L}: m0={s['m0']} xt={s['x_t']:.3f} xs={s['x_s']:.3f} xc={s['x_c']:.3f} g1={s['gap1']:.3f} |"
        else:
            line += f" L{L}: m0={s['m0']} |"
    if len(Ls) >= 2 and 'x_c' in S[Ls[-1]]:
        L1, L2 = Ls[-2], Ls[-1]
        v = S[L2]['x_c'] / (2 * np.pi)
        cv = 6 * (S[L2]['E0'] / L2 - S[L1]['E0'] / L1) / (np.pi * (1 / L1**2 - 1 / L2**2))
        line += f" r_s={S[L2]['x_s']/S[L2]['x_t']:.3f} r_c={S[L2]['x_c']/S[L2]['x_t']:.3f} v={v:.3f} c_cas={cv/v if v>1e-9 else np.nan:.3f}"
        line += f" low{L2}={S[L2]['low'][:4]}"
    print(line)

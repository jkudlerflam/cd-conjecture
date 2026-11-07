"""Probe refined landing points: spectra at L = 7..14 and a one-parameter cut that rescales the field,
x(s) = (J, (1+s) B). A band-onset (PT) point shows a z = 2 gap at s = 0, a gap on one side and a
c = 1 liquid (L*gap of order one, many low states) on the other."""
import json, sys
import numpy as np
import search
from ed import bond, spectrum, all_levels, params_to_JB, gauge_M
fn = sys.argv[1]; idx = int(sys.argv[2])
d = [json.loads(l) for l in open(fn)][idx]
x = np.array(d["x"]); J, B = params_to_JB(x)
M = gauge_M(bond(J, B))[0]
print(f"point {idx}: g={d['g_final']:.3f}, M={M:.4f}, J=\n{np.round(J,4)}\nB={np.round(B,4)}")
for L in [7, 8, 9, 10, 11, 12, 13, 14]:
    w, k = all_levels(spectrum(bond(J, B), L)); o = np.argsort(w); w = w[o]; k = k[o]
    e = w - w[0]
    print(f"  L={L:2d}: n(<M/L)={int(np.sum(e < M/L)):3d}  L*e1..e4={np.round(L*e[1:5],3).tolist()}  L^2*e1={L*L*e[1]:.3f}  k of lowest 6: {k[:6].tolist()}", flush=True)
print("cut x(s) = (J, (1+s)B), L = 10, 12:")
for s in [-0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2]:
    row = []
    for L in [10, 12]:
        w, _ = all_levels(spectrum(bond(J, (1 + s) * B), L)); w = np.sort(w); e = w - w[0]
        row.append((int(np.sum(e < M / L)), round(L * e[1], 3)))
    print(f"  s={s:+.2f}: (n(<M/L), L*e1) at L=10,12: {row}", flush=True)

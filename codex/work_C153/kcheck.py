#!/usr/bin/env python3
"""Compare ring ED K, v with the near-boundary prediction K0 = (pi/2)/sqrt(r(1+r)(r-u)), v0 = 4 sqrt(r(r-u)/(1+r))."""
import numpy as np, sys
from ed import luttinger
Ls = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [12, 16, 20]
for b in [0.5, 1.0, 2.0, 4.0, 8.0]:
    r = np.sqrt(1 + b * b / 4)
    for w in [0.002, 0.01, 0.05, 0.2]:
        u = r - w
        K0 = (np.pi / 2) / np.sqrt(r * (1 + r) * w); v0 = 4 * np.sqrt(r * w / (1 + r))
        out = [luttinger(L, u, b) for L in Ls]
        print(f'b={b} r-u={w} K0={K0:.3f} v0={v0:.4f} | ' + ' '.join(f'L{d["L"]}:K={d["K"]:.3f},v={d["v"]:.4f}' for d in out), flush=True)

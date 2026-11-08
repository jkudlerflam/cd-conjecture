#!/usr/bin/env python3
"""Ring ED scan of the Luttinger parameter K(L) (charge gap and twist) along cuts of fixed b, L = 12, 16, 20.
Where the chain is gapped, K(L) falls with L toward 0; in the liquid it is L-independent and > 2.
The crossing K(L) = 2 with weak L dependence gives a rough BKT estimate (log corrections not removed).
Written for codex/work_C153 (Claude Opus agent, 2026-11-07)."""
import numpy as np, json, sys
from ed import luttinger
res = []
for b in [0.25, 0.5, 1.0, 1.5, 2.0, 3.0]:
    r = np.sqrt(1 + b * b / 4)
    us = sorted(set([round(x, 4) for x in list(np.linspace(-0.8, 0.6, 8)) + list(r - np.array([0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.01]))]))
    for u in us:
        if u <= -1 or u >= r: continue
        row = [luttinger(L, u, b) for L in (12, 16, 20)]
        res.append(dict(b=b, u=u, r=r, K=[d['K'] for d in row], v=[d['v'] for d in row], gap1=[d['dE1'] * d['L'] for d in row]))
        print(f'b={b} u={u:.4f} r-u={r-u:.4f} K(12,16,20)=' + ','.join(f'{d["K"]:.3f}' for d in row) + '  L*dE1=' + ','.join(f'{d["dE1"]*d["L"]:.3f}' for d in row), flush=True)
json.dump(res, open('kscan.json', 'w'))

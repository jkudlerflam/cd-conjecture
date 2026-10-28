"""Resonance margins of the bare second-order bands, with the minimising
four-fermion process identified (n = 2 at t_d = 1e-3, n = 3 at t_d = 1e-5 and 0)."""
import json, itertools, numpy as np
from rt_bond import P
from rt_cell import second_order, hs_quadratic, herm_part_analysis, bands
def d(x):
    return abs((x + np.pi) % (2 * np.pi) - np.pi)
out = []
for n, td in ((2, 1e-3), (2, 0.0), (3, 1e-5), (3, 1e-4), (3, 0.0)):
    H2 = second_order(n, P(td=td, g=0.05))['H2']
    e0, h, res, _ = hs_quadratic(H2, n)
    T, mu = h[:n, n:], h[:n, :n] + h[n:, n:]
    ev, cr, asym = herm_part_analysis(T)
    sc = np.abs(T).max()
    eps = 0.5 * ev[0] * sc
    bs = bands(T, mu + eps * np.eye(n), nk=400001)
    pts = [(b, r) for b, B in enumerate(bs) for r in B['roots']]
    best = (np.inf, None)
    pairs = list(itertools.combinations_with_replacement(range(len(pts)), 2))
    for (a, b) in pairs:
        for (c, e) in pairs:
            if (a, b) == (c, e):
                continue
            # skip trivially forward combinations {a,b} = {c,e} handled above; also skip same-band same-point
            m = d(pts[a][1] + pts[b][1] - pts[c][1] - pts[e][1])
            if m < best[0]:
                best = (m, ((pts[a][0], round(pts[a][1], 4)), (pts[b][0], round(pts[b][1], 4)),
                            (pts[c][0], round(pts[c][1], 4)), (pts[e][0], round(pts[e][1], 4))))
    r = dict(n=n, td=td, centres=[B['centre'] for B in bs], fillings=[B['filling'] for B in bs],
             n_roots=[len(B['roots']) for B in bs], delta4=best[0], argmin=best[1])
    print(json.dumps(r, default=str), flush=True); out.append(r)
json.dump(out, open('margins.json', 'w'), indent=1, default=str)

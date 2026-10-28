"""Largest loop link t_d keeping the Hermitian part of the rotated T positive,
against D_x (leading order in g).  The loop gives the right port direct,
unsuppressed access to the far end of the stub, so t_d* should fall with the
express suppression."""
import json, numpy as np
from rt_bond import P
from rt_token import T_token, crawford
out = []
for n in (2, 3, 4):
    for Dx in (4.0, 10.0, 30.0):
        if n == 4 and Dx == 4.0:
            continue
        tds = np.logspace(-12, 0, 49)
        last = None
        for td in tds:
            T = T_token(n, P(td=td, Dx=Dx, DM=Dx))
            cr, th = crawford(T, nth=1441)
            if cr > 0:
                last = td
            else:
                break
        T0 = T_token(n, P(td=0.0, Dx=Dx, DM=Dx))
        r = dict(n=n, Dx=Dx, td_star=last, T0_scale=float(np.abs(T0).max()))
        out.append(r); print(r, flush=True)
json.dump(out, open('scan_td.json', 'w'), indent=1)

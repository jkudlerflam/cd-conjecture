"""Delta > 0 (needed so that long cells keep the code as ground space): does the
Hermitian part of the rotated T stay positive?  Leading order in g, t_d = 0,
neutral express phase pi/(2n)."""
import json, numpy as np
from rt_bond import P
from rt_token import T_token, rs_eigs
out = []
for n in range(2, 9):
    for Delta in (0.0, 0.1, 0.3, 1.0):
        for Dx in (10.0, 30.0, 100.0):
            p = P(td=0.0, Delta=Delta, Dx=Dx, DM=Dx, ph_x={'0': np.pi / (2 * n), 'p': 0.0, 'm': 0.0})
            T = T_token(n, p)
            ev = rs_eigs(T)
            r = dict(n=n, Delta=Delta, Dx=Dx, T_scale=float(np.abs(T).max()), lam_min=float(ev[0]))
            out.append(r)
        print(n, Delta, [('%g' % r['Dx'], '%.2e' % r['lam_min'], '%.1e' % r['T_scale']) for r in out[-3:]], flush=True)
json.dump(out, open('scan_delta.json', 'w'), indent=1)

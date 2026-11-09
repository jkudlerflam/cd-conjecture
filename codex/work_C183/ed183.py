#!/usr/bin/env python3
"""C183: fc.py (--cp path, dmrg2.py) against sparse ED (c166lib.ed_e0) at L = 10, 12 for bonds C183 adds:
off the stratum (J_x + J_y != 0, complex path), with XZ - ZX (Kxz), with XZ + ZX and h_x, at r = 4 and 6.  Also the
parity-sector option against ED restricted by parity (energy of the lowest odd state)."""
import json, subprocess, numpy as np
import c166lib as X, tune183 as TU
cases = []
for r, wa, phi, extra in ((4.0, 2.0, np.pi / 2, {}), (6.0, 1.0, 1.50, {}), (6.0, 1.0, np.pi / 2, dict(Kxz=0.3)),
                          (6.0, 2.0, np.pi / 2, dict(Jxz=0.3, hx=0.3))):
    T = TU.tune(r, wa, phi)
    ex = {k: v * T['t2'] for k, v in extra.items()}
    cp = TU.couplings(T, 0.1, ex)
    cases.append((f"r{r:g}w{wa:g}p{phi:.2f}" + ''.join(extra), cp))
for tag, cp in cases:
    jxz, kxz, hx = cp.get('Jxz', 0.0), cp.get('Kxz', 0.0), cp.get('hx', 0.0)
    J = np.array([[cp['Jx'], cp['D'], jxz + kxz], [-cp['D'], cp['Jy'], 0], [jxz - kxz, 0, cp['Jz']]]); f = np.array([hx, 0, cp['hz']])
    for L in (10, 12):
        e_ed = X.ed_e0(J, f, L)[0]
        args = ['python3', 'fc.py', f'chk_{tag}', str(L), '--cp', json.dumps(cp), '--chis', '16', '32', '64', '64', '64', '--outdir', 'edchk', '--nocorr']
        subprocess.run(args, capture_output=True)
        d = json.load(open(f'edchk/chk_{tag}_{L}.json'))
        print(tag, L, f'ED {e_ed:.12f} DMRG {d["E"]:.12f} diff {d["E"]-e_ed:.1e} gauge={d["gauge"]}', flush=True)

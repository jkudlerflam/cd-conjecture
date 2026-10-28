"""Gram-limit test: lambda_min of the Hermitian part of the phase-rotated T
(leading order in g, t_d = 0) against D_x, and the bare-stub Gram matrix."""
import json, numpy as np
from rt_bond import P
from rt_token import T_token, rs_eigs, crawford
out = []
for n in range(2, 8):
    G = T_token(n, P(td=0.0), stub_only=True).real
    gev = np.linalg.eigvalsh((G + G.T) / 2) / np.abs(G).max()
    rec = dict(n=n, gram_eigs=list(gev), gram_asym=float(np.abs(G - G.T).max() / np.abs(G).max()), rows=[])
    print('n', n, 'Gram eig min/max %.3e' % (gev[0] / gev[-1]), 'asym', rec['gram_asym'], flush=True)
    for Dx in (2.0, 4.0, 10.0, 30.0, 100.0, 1000.0):
        p = P(td=0.0, Dx=Dx, DM=Dx)
        T = T_token(n, p)
        ev = rs_eigs(T)
        cr, th = crawford(T, nth=3601)
        R = T * np.exp(-1j * np.angle(T[0, 0]))
        Rn = (R / np.abs(R).max()).real
        dist = float(np.abs(Rn / np.sign(Rn[0, 0]) - G / np.abs(G).max()).max())
        row = dict(Dx=Dx, T_scale=float(np.abs(T).max()), lam_min=float(ev[0]), lam_min_over_max=float(ev[0] / ev[-1]),
                   crawford=float(cr), imag_after_phase=float(np.abs(R.imag).max() / np.abs(R).max()), dist_to_gram=dist)
        rec['rows'].append(row)
        print('   Dx %6g |T| %.3e  lam_min/|T| %+.3e  crawford %+.3e  Im/|T| %.1e  |T/|T| - Gram| %.2e' %
              (Dx, row['T_scale'], row['lam_min'], cr, row['imag_after_phase'], dist), flush=True)
    out.append(rec)
json.dump(out, open('scan_gram.json', 'w'), indent=1)

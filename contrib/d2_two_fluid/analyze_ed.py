#!/usr/bin/env python3
"""Summary of the exact ring check of TDA survivors (data/ed_survivors_L*.json)."""
import json, sys, numpy as np
for L in (12, 14):
    R = json.load(open(f'data/ed_survivors_L{L}.json'))
    pre = dense = dil = 0; rows = []
    for r in R:
        e = r['ed']; p = r['p']
        if 'mu_ex' not in e:
            continue
        if e.get('below') is None:
            pre += 1; rows.append((p, r['mu_star_tda'], e['mu_ex'], e['N_gs'], None)); continue
        if e['N_gs'] / L > 0.3:
            dense += 1
        else:
            dil += 1
        exc = np.array(e['below']['exc']); E, N = exc[:, 0], exc[:, 1]
        vals = [k for k in range(L) if E[k] <= E[k - 1] and E[k] <= E[(k + 1) % L]]
        vs = sorted([(E[k], N[k], k) for k in vals])
        rows.append((p, r['mu_star_tda'], e['mu_ex'], e['N_gs'], vs[:3]))
    print(f"L={L}: {len(R)} survivors; ground state already non-vacuum at mu*_tda-0.6: {pre}; "
          f"onset into a dense state (<N>/L > 0.3): {dense}; dilute onset: {dil}")
    if L == 14 or '-v' in sys.argv:
        for p, ms, mex, Ngs, vs in rows:
            s = ' '.join(f"(E={a:.3f},N={b:.2f},k={c})" for a, b, c in vs) if vs else 'preempted'
            print(f" phi={p['phi']:.2f} V={p['V']:+.2f} gL={p['gL']:.2f} |gR|={np.hypot(*p['gR']):.2f} "
                  f"mu*_tda={ms:.2f} mu_ex={mex:.2f} N_gs={Ngs:.2f} {s}")

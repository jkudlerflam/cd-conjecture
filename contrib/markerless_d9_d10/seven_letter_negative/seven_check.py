"""Negative result: the 7-letter reflection design (no express rail; both tokens enter from the port ends).
Its hopping is T = (reflection) x (positive Gram): dimerised, at most n mod 2 bands cross zero."""
import numpy as np, sys, json
from twocell import analyse, one_particle, jw_many_body
from bands import fermi_roots
out = []
for n in (2, 3):
    for label, gL in [('one rung', 0.0), ('two rungs with flux', 0.6*np.exp(1.1j))]:
        p = dict(lam=1.0, gR=1.0, gL=gL, DA=0.7, DB=0.7, nu=0.0, f=0.0, chi=0.0)
        r = analyse(n, p); h, vacE = one_particle(r['h2'], n, r['bitlist'])
        M = jw_many_body(h, vacE, r['bitlist'])
        T = h[:n,n:]; mu = h[:n,:n] + h[n:,n:]
        z = fermi_roots(mu, T); on = int((np.abs(np.abs(z)-1) < 1e-8).sum())
        J = np.eye(n)[::-1]; S = J@T
        rec = dict(n=n, design=label, leaks=r['leaks'], gaussian_residual=float(np.abs(r['h2']-M).max()/np.abs(T).max()),
                   mu_over_T=float(np.abs(mu).max()/np.abs(T).max()), fermi_points=on, needed=2*n,
                   moduli=np.round(np.sort(np.abs(z)),4).tolist())
        out.append(rec); print(rec, flush=True)
json.dump(out, open('seven_check.json','w'), indent=1)

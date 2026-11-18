"""Exact two-cell checks of the 9-letter bond at n = 2, 3 (well-conditioned diagnostic parameters) and
agreement of the arbitrary-precision scalar formula with the brute-force letter model."""
import numpy as np, json, mpmath as mp
from twocell9 import analyse, one_particle, jw_many_body
from scalar9 import blocks
mp.mp.dps = 40
out = []
for n in (2, 3):
    for a, D, td in [(2.6, 3.0, 0.07), (2.2, 0.4, 0.3), (3.0, 8.0, 0.0)]:
        p = dict(lam=1.0, g=1.0, Delta=a-2*np.cos(np.pi/(n+2)), D=D, nu=0.0, f=0.0, phiA=np.pi/(2*(n+1)), td=td)
        r = analyse(n, p); h, vacE = one_particle(r['h2'], n, r['bitlist'])
        M = jw_many_body(h, vacE, r['bitlist'])
        al4 = (2/(n+2))**2*np.sin(np.pi/(n+2))**4
        T, muL, muR = blocks(n, a=mp.mpf(a), D=mp.mpf(D), td=mp.mpf(td))
        c = lambda X: np.array(X.tolist(), dtype=complex)*al4
        sc = np.abs(h[:n,n:]).max(); h0 = r['h0']
        rec = dict(n=n, a=a, D=D, td=td, two_cell_states=r['nstates'], leaks=r['leaks'], code_residual=float(r['code_resid']),
            first_order=float(r['first']), g0_bit_dependence=float(np.ptp(np.diag(h0).real)), g0_offdiagonal=float(np.abs(h0-np.diag(np.diag(h0))).max()),
            g1_coefficient=float(np.abs(r['h1']).max()), g2_hermiticity=float(np.abs(r['h2']-r['h2'].conj().T).max()/sc),
            gaussian_residual_over_T=float(np.abs(r['h2']-M).max()/sc),
            scalar_vs_bruteforce_T=float(np.abs(h[:n,n:]-c(T)).max()/sc), scalar_vs_bruteforce_muL=float(np.abs(h[:n,:n]-c(muL)).max()/sc),
            scalar_vs_bruteforce_muR=float(np.abs(h[n:,n:]-c(muR)).max()/sc), mu_over_T=float(np.abs(h[:n,:n]+h[n:,n:]).max()/sc), maxT=float(sc))
        out.append(rec); print(rec, flush=True)
json.dump(out, open('two_cell_checks.json','w'), indent=1)

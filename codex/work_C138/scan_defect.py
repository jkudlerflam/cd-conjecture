"""C138 Task 5 scan: doubling defect and bulk deviation near Delta=-1.

For eta with epsilon = 1+Delta = 10^-k, k=2..6, rings L=2^m (m=2..13) and
sectors n = nu L (dyadic nu) plus a few dilute and near-half sectors.
Writes energies G_L(n) = F_L(n) - Delta L/4 to energies_<k>.json (cache),
and derived ratios to scan_<k>.json. Floating point only.

Usage: python3 scan_defect.py K [LMAX]
"""
import json
import sys
import time
from pathlib import Path
import numpy as np
from bethe_ferro import solve, Model

K = int(sys.argv[1])
LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 8192
eps = 10.0**(-K)
eta = 2*np.arcsin(np.sqrt(eps/2))
m = Model(eta)
cache_file = Path(f'energies_{K}.json')
cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
NUS = [1/2, 7/16, 3/8, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128]


def G(L, n):
    key = f'{L},{n}'
    if key not in cache:
        t = time.time()
        g, x, d = solve(L, n, eta)
        cache[key] = dict(G=g, res=d['res'], steps=d['steps'],
                          xmax=float(x[-1]) if n else 0.0,
                          mingap=float(np.min(np.diff(x))) if n > 1 else None)
        if time.time()-t > 5:
            cache_file.write_text(json.dumps(cache))
    return cache[key]['G']


Ls = [2**j for j in range(2, 30) if 2**j <= LMAX//2]
rows = []
for L in Ls:
    ns = set([1, 2, 3, 4, 8, int(np.sqrt(L)), L//2, L//2-1, L//2-2,
              L//2-int(np.sqrt(L))])
    for nu in NUS:
        if nu*L == int(nu*L):
            ns.add(int(nu*L))
    for n in sorted(ns):
        if not 1 <= n <= L//2:
            continue
        g1, g2 = G(L, n), G(2*L, 2*n)
        lam = L*(g2-2*g1)       # L[F_2L(2n)-2F_L(n)] in units of J
        Q = 2*lam               # in units of M_0 = J/2
        rows.append(dict(L=L, n=n, nu=n/L, Lambda_d=lam, Q=Q,
                         band=eps*L*L, dilute=2*np.pi**2*n**3/L,
                         Q_over_eta=Q/eta, etaL=eta*L, epsL2=eps*L*L))
    cache_file.write_text(json.dumps(cache))
    print(f'K={K} L={L} done', flush=True)

# Ladder estimate of L(F_L - L e_bulk) = -(1/2) sum_j Lambda_j / 4^j.
lam_of = {(r['L'], r['n']): r['Lambda_d'] for r in rows}
ladder = []
for (L, n) in sorted(lam_of):
    terms = []
    LL, nn = L, n
    while (LL, nn) in lam_of:
        terms.append(lam_of[(LL, nn)])
        LL, nn = 2*LL, 2*nn
    if len(terms) < 3:
        continue
    s = -0.5*sum(t/4**j for j, t in enumerate(terms))
    tail = 0.5*max(abs(t) for t in terms)*4.0**(-len(terms))*4/3
    ladder.append(dict(L=L, n=n, nu=n/L, dev=s, tail_bound_if_bounded=tail,
                       terms=len(terms), dev_over_eta=s/eta))

v = (np.pi/2)*np.sin(eta)/(np.pi-eta)
out = dict(K=K, eps=eps, eta=eta, v=v, casimir_Lambda=np.pi*v/4,
           casimir_dev=-np.pi*v/6, rows=rows, ladder=ladder,
           max_Q=max(r['Q'] for r in rows),
           max_absQ=max(abs(r['Q']) for r in rows),
           min_Q=min(r['Q'] for r in rows),
           max_abs_dev=max(abs(r['dev']) for r in ladder) if ladder else None,
           band_violations=sum(1 for r in rows if r['Q'] > r['band']+1e-12),
           dilute_violations=sum(1 for r in rows if r['Q'] > r['dilute']+1e-12))
Path(f'scan_{K}.json').write_text(json.dumps(out, indent=1))
print(json.dumps({k: out[k] for k in ['K', 'eta', 'casimir_Lambda', 'max_Q', 'min_Q',
                                      'max_abs_dev', 'band_violations',
                                      'dilute_violations']}))

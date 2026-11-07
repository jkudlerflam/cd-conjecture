"""Compact tables from analysis_majorant: V'_min/M for S=2 and kappa_eff, restricted to Trotter step beta/m <= dmax."""
import numpy as np, sys
from analysis_majorant import summary, MODELS
dmax = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
for Delta in MODELS:
    print('Delta=%g M=%g  (beta/m <= %g)' % (Delta, MODELS[Delta], dmax))
    print('   m | V/M S=2,K=ln2 (L,beta) | S=2,K=2 | S=6,K=ln2 | kappa(K=0) kappa(K=ln2) | max logZt')
    for r in summary(Delta, range(2, 10), Ss=(2, 6), Ks=(np.log(2), 2.0), dmax=dmax):
        a = r['S2_K0.69']; b = r['S2_K2.00']; c = r['S6_K0.69']
        print('  %2d | %7.3f (%d,%g) | %7.3f (%d,%g) | %7.3f | %7.3f %7.3f | %6.3f' % (r['m'], a[0], a[1], a[2], b[0], b[1], b[2], c[0], r['kappa_K0.00'], r['kappa_K0.69'], r['maxlogZt']))

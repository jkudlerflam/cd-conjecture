"""C161 (Claude Opus agent, 2026-11-08). Classical ground state of the K-stratum functional of soliton_valley.py:
planar CSL vs tilted CSL (one soliton per cell, period optimized) vs uniform +x (energy 0).
"""
import sys
import numpy as np
from soliton_pair import ring_min

kappa, mu = 0.05, 2.3
print(f'kappa={kappa} mu={mu}')
print('h       best_planar(e,T)          best_tilted(e,T,max|n_z|, <n_z>)   ground state')
for h in [float(a) for a in sys.argv[1:]] or (0.0, 0.01, 0.02, 0.03):
    bp = bt = None
    for T in (5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0):
        ep, npl = ring_min(h, kappa, mu, T, 1, [0.0], dx=0.1)
        et, nt = ring_min(h, kappa, mu, T, 1, [0.6], dx=0.1)
        if bp is None or ep / T < bp[0]:
            bp = (ep / T, T)
        if bt is None or et / T < bt[0]:
            bt = (et / T, T, np.abs(nt[:, 2]).max(), nt[:, 2].mean())
    gs = 'uniform' if min(bp[0], bt[0]) > -1e-7 else ('tilted CSL' if bt[0] < bp[0] - 1e-7 and bt[2] > 0.05 else 'planar CSL')
    print(f'{h:.3f}  ({bp[0]:+.5f}, {bp[1]:4.1f})   ({bt[0]:+.5f}, {bt[1]:4.1f}, {bt[2]:.3f}, {bt[3]:+.3f})   {gs}')
    sys.stdout.flush()

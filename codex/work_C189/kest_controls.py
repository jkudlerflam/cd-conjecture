"""C189 Part 1 calibration of the K estimators on the controls (C176's c1 and c2 runs, exact ff2v open chains).
Usage: python3 kest_controls.py > kest_controls.out"""
import numpy as np, warnings
from kest import *
warnings.simplefilter('ignore')
D = '../work_C176/typhon_out/out/'
Kc1 = (1 - 0.25) ** 2
print('E_c: K from S(q) on (0.05,0.3) [neighbours (0.05,0.2), (0.1,0.3)]. Exact: c1 0.5625 (spinless), c2 1, ff2v 1')
for L in (64, 128, 192, 256):
    c1 = load_run(D + f'c1_L{L}.jsonl'); c2 = load_run(D + f'c2_L{L}.jsonl'); ff = ff_run(L, 3 * L // 16)
    s = f'  L={L}:'
    for nm, run, nc in (('c1', c1, 1), ('c2', c2, 2), ('ff2v', ff, 2)):
        v = [K_from_Sq(run, lo, hi, nc)[0] for lo, hi in ((0.05, 0.3), (0.05, 0.2), (0.1, 0.3))]
        s += f"  {nm} {v[0]:.3f} [{v[1]:.3f}, {v[2]:.3f}]"
    print(s)
print("E_F': Friedel exponent on (8, L/4) [neighbours (8, L/2), (12, L/4)]. Exact: c1 0.5625 at 2kF=2 pi n; c2 1 at pi/4, 3pi/4; ff2v 1 at pi n (+ 2k0 term)")
for L in (64, 128, 192, 256, 512):
    ff = ff_run(L, 3 * L // 16)
    rows = [('ff2v', ff, [np.pi * n_mid(ff), 2 * K0])]
    if L <= 256:
        c1 = load_run(D + f'c1_L{L}.jsonl'); c2 = load_run(D + f'c2_L{L}.jsonl')
        rows = [('c1', c1, [2 * np.pi * n_mid(c1)]), ('c2', c2, [np.pi / 4, 3 * np.pi / 4])] + rows
    s = f'  L={L}:'
    for nm, run, ks in rows:
        v = [friedel_fit(run, lo, hi, ks)['eta'] for lo, hi in ((8, L // 4), (8, L // 2), (12, L // 4))]
        s += f"  {nm} {v[0]:.3f} [{v[1]:.3f}, {v[2]:.3f}]"
    print(s)
print('E_D: pair envelope exponent on (16, L/4) [neighbours (24, L/4), (16, 3L/8)], chord distance. Exact ff2v 2')
for L in (128, 192, 256, 512):
    ff = ff_run(L, 3 * L // 16)
    v = [eta(ff, 'D', lo, hi) for lo, hi in ((16, L // 4), (24, L // 4), (16, 3 * L // 8))]
    print(f'  L={L}: ff2v {v[0]:.3f} [{v[1]:.3f}, {v[2]:.3f}]')
print('E_G (excluded, fails controls): one-body envelope exponent on (16, 64), chord. Exact c1 1.170, c2 1.5, ff2v 1')
for L in (128, 192, 256):
    c1 = load_run(D + f'c1_L{L}.jsonl'); c2 = load_run(D + f'c2_L{L}.jsonl'); ff = ff_run(L, 3 * L // 16)
    print(f'  L={L}: c1 {eta(c1, "G", 16, 64):.3f}  c2 {eta(c2, "G", 16, 64):.3f}  ff2v {eta(ff, "G", 16, 64):.3f}')
print('Friedel spectrum of ff2v at L = 1024 (oscillating part weighted by the chord factor): top peaks (k, |FT|)')
L = 1024; ff = ff_run(L, 3 * L // 16); dn = ff['dens']; nb = n_mid(ff); x = np.arange(L)
osc = (dn - nb) * np.sin(np.pi * (x + 1) / (L + 1))
ks = np.linspace(0, np.pi, 4001); F = np.abs(np.exp(1j * np.outer(ks, x[:L // 2])) @ osc[:L // 2])
idx = sorted([i for i in range(1, len(ks) - 1) if F[i] > F[i - 1] and F[i] > F[i + 1]], key=lambda i: -F[i])[:4]
print('  pi n =', round(np.pi * nb, 4), ' 2k0 folded =', round(2 * np.pi - 2 * K0, 4), ' 2(k0+-kF) folded =',
      [round(min(k % (2 * np.pi), 2 * np.pi - k % (2 * np.pi)), 4) for k in ks_two_valley(nb)], ' peaks:',
      [(round(ks[i], 4), round(F[i], 3)) for i in idx])

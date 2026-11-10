"""C189 Part 1, boson-frame estimators (amendment A3, post hoc, after the coordinator's frame check).
Relations (PROOF 3): lattice boson b at a valley momentum ~ e^{i theta_a}: exponent eta_B = 1/(4K_c) + 1/(4K_s);
momentum-pi inter-valley pair b_j b_{j+d} (d = 2, 3; equal to -c_j c_{j+d} under the blockade) ~ e^{i sqrt2 theta_c}:
exponent eta_P = 1/K_c. Hence K_c^P = 1/eta_P and 1/K_s^B = 4 eta_B - eta_P.
Fits: rows from a fixed site with the open-chain chord distance and the image factor (PROOF 5), centred RMS
envelope over +-4 samples. Windows: boson (8, L/2-12) [neighbours (16, L/2-12), (8, L/4)], rows starting 4 to 10
sites inside the central window; pair (16, 3L/8) [neighbours (24, 3L/8), (16, L/4)].
Controls: c2 (boson 0.5, inter-leg pair d = 3: 1, intra-leg pair d = 2: 2), c1 (spinless, boson 1/(2K) = 0.889).
Usage: python3 kest_boson.py > kest_boson.out"""
import numpy as np, warnings
from kest import *
warnings.simplefilter('ignore')
D = '../work_C176/typhon_out/out/'


def gb_rows(run, lo, hi, starts=(4, 6, 8, 10)):
    G = run['Gb']; a0 = run['gwin'][0]; L = run['L']; v = []
    for s in starts:
        row = np.abs(G[s, s + 1:]); r = np.arange(1, 1 + len(row))
        v.append(fit_chord(r, row, L, lo, hi, 'row', a0 + s, image=True))
    return float(np.mean(v))


def pairfit(run, sep, lo, hi):
    L = run['L']; v = []
    for i0, row in D_rows(run, sep).items():
        r = np.arange(sep + 1, sep + 1 + len(row)); v.append(fit_chord(r, row, L, lo, hi, 'row', i0, image=True))
    return float(np.mean(v))


def est(run):
    L = run['L']
    b = [gb_rows(run, lo, hi) for lo, hi in ((8, L // 2 - 12), (16, L // 2 - 12), (8, L // 4))]
    p2 = [pairfit(run, 2, lo, hi) for lo, hi in ((16, 3 * L // 8), (24, 3 * L // 8), (16, L // 4))]
    p3 = [pairfit(run, 3, lo, hi) for lo, hi in ((16, 3 * L // 8), (24, 3 * L // 8), (16, L // 4))]
    return b, p2, p3


print('Controls')
bias = {}
for L in (128, 192, 256):
    c2 = load_run(D + f'c2_L{L}.jsonl'); c1 = load_run(D + f'c1_L{L}.jsonl')
    b2, p22, p23 = est(c2); b1, _, _ = est(c1)
    bias[L] = dict(B=max(abs(b2[0] - 0.5), abs(b1[0] - 8 / 9)), P=abs(p23[0] - 1.0))
    print(f'  L={L}: c2 eta_B={b2[0]:.3f} (0.5)  pair d=3 {p23[0]:.3f} (1)  pair d=2 {p22[0]:.3f} (2);  c1 eta_B={b1[0]:.3f} (0.889)'
          f'  -> control bias B {bias[L]["B"]:.3f}, P {bias[L]["P"]:.3f}')
print('\nCandidates (phi = pi/2, |w| = 1): eta_B, eta_P (d = 2; d = 3 in brackets), K_c^P = 1/eta_P, K_s^B = 1/(4 eta_B - eta_P)')
for s in ('sym316', 'sym532', 'sym18'):
    for L in (128, 192, 256):
        run = load_run(D + f'{s}_L{L}.jsonl')
        b, p2, p3 = est(run)
        sB = np.hypot(0.5 * (max(b) - min(b)), bias[L]['B'])
        sP = np.hypot(0.5 * (max(p2) - min(p2)), bias[L]['P'])
        inv = 4 * b[0] - p2[0]; sinv = np.hypot(4 * sB, sP)
        Ks = 1 / inv
        lo2, hi2 = 1 / (inv + 2 * sinv), (1 / (inv - 2 * sinv) if inv - 2 * sinv > 0 else np.inf)
        print(f'  {s} L={L}: eta_B={b[0]:.3f}+-{sB:.3f}  eta_P={p2[0]:.3f}+-{sP:.3f} [{p3[0]:.3f}]  K_c^P={1 / p2[0]:.3f}+-{sP / p2[0] ** 2:.3f}'
              f'  1/K_s={inv:.3f}+-{sinv:.3f}  K_s^B={Ks:.3f}  (2 sigma range {lo2:.2f} to {hi2:.2f})')

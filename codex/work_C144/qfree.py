"""Doubling ratio Q_L = L (E_2L - 2 E_L) for a gapped (dimerized) Majorana ring, hoppings 1 +- d, best BC sector,
compared with the critical value 2 L tan(pi/4L)/cos(pi/2L) (XZ, closed form, PROOF 2.1). Checks the
'no Casimir overshoot' input of the Q1 argument for a free massive sector."""
import numpy as np
def E(L, d):
    best = 0
    for bc in (-1, 1):
        q = 2 * np.pi * (np.arange(L // 2) + (0.5 if bc == -1 else 0.0)) / (L // 2)
        e = -np.sum(np.abs((1 + d) + (1 - d) * np.exp(1j * q)))
        best = min(best, e)
    return best
for L in (4, 8):
    print("check XZ closed form L=%d:" % L, L * (E(2 * L, 0) - 2 * E(L, 0)), 2 * L * np.tan(np.pi / (4 * L)) / np.cos(np.pi / (2 * L)))
for d in (0.0, 0.01, 0.05, 0.2):
    Q = [L * (E(2 * L, d) - 2 * E(L, d)) for L in (4, 8, 16, 32, 64, 128, 256, 512)]
    print(f"d={d}: Q_L (L=4..512) =", np.round(Q, 4).tolist(), " max =", round(max(Q), 4))

"""C186: check ring_mixed.py (momentum ED) against a direct real-space construction of the same model, L = 12, 13."""
import numpy as np
from ring_mixed import Ring, model, pconfigs
def realspace(L, m):
    S = pconfigs(L); pos = {int(s): i for i, s in enumerate(S)}
    H = np.zeros((len(S), len(S)), complex)
    o = lambda s, j: (int(s) >> (j % L)) & 1
    for i, s in enumerate(S):
        s = int(s)
        for j in range(L):
            if o(s, j):
                H[i, i] += -m['mu'] + m['V2'] * o(s, j + 2) + m['V3'] * o(s, j + 3)
                b = s & ~(1 << j)
                for dj, a in ((1, m['w']), (-1, np.conj(m['w']))):
                    if not o(s, j + 2 * dj): H[pos[b | (1 << ((j + dj) % L))], i] += a
                for dj in (2, -2):
                    sg = 1 if dj > 0 else -1
                    if not o(s, j + dj) and not o(s, j + dj + sg): H[pos[b | (1 << ((j + dj) % L))], i] += -1
                amp = m['lam'] + m['r'] * (o(s, j + 2) - o(s, j - 2))
                H[pos[b], i] += amp
            elif not o(s, j - 1) and not o(s, j + 1):
                amp = m['lam'] + m['r'] * (o(s, j + 2) - o(s, j - 2))
                H[pos[s | (1 << j)], i] += amp
    return H
for L in (12, 13):
    m = model(1.3, 1.45, 0.4, -2.3, lam=0.07)
    H = realspace(L, m)
    assert abs(H - H.conj().T).max() < 1e-12
    e_rs = np.sort(np.linalg.eigvalsh(H))
    ring = Ring(L)
    e_k = []
    for kk in range(L):
        Hk, N = ring.H(m, 2 * np.pi * kk / L)
        e_k += list(np.linalg.eigvalsh(Hk.toarray()))
    e_k = np.sort(e_k)
    print(f"L={L}: dim {len(e_rs)} vs {len(e_k)}, max |diff| = {np.abs(e_rs - e_k).max():.2e}")

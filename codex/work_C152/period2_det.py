"""Check of Lemma P2 (PROOF section 3): for a period-2 nearest-neighbour BdG chain, det H_BdG(k) is a polynomial of
degree <= 2 in cos k (Fourier modes |n| <= 2 only, and even in k). Prints the largest |Fourier coefficient| with
|n| >= 3 and the largest odd (sin) component, relative to the largest coefficient, over random models."""
import numpy as np
from period2_mf import model, K
rng = np.random.default_rng(3)
def detk(p):
    hR, DR = model(p)
    h = sum(m[None] * np.exp(1j * K * R)[:, None, None] for R, m in hR.items())
    hm = sum(m[None] * np.exp(-1j * K * R)[:, None, None] for R, m in hR.items())
    D = sum(m[None] * np.exp(1j * K * R)[:, None, None] for R, m in DR.items())
    top = np.concatenate([h, D], axis=2)
    bot = np.concatenate([np.conj(np.transpose(D, (0, 2, 1))), -np.transpose(hm, (0, 2, 1))], axis=2)
    return np.linalg.det(np.concatenate([top, bot], axis=1))
worst_hi = 0; worst_odd = 0; worst_im = 0
for i in range(500):
    p = rng.normal(size=10); d = detk(p)
    worst_im = max(worst_im, np.max(np.abs(d.imag)) / np.max(np.abs(d)))
    F = np.fft.fft(d.real) / len(K)
    n = np.fft.fftfreq(len(K), 1 / len(K))
    big = np.max(np.abs(F))
    worst_hi = max(worst_hi, np.max(np.abs(F[np.abs(n) >= 3])) / big)
    worst_odd = max(worst_odd, np.max(np.abs(F - F[(-np.arange(len(K))) % len(K)])) / big)
print('max rel |Im det| =', worst_im, ' max rel Fourier |n|>=3 =', worst_hi, ' max rel odd part =', worst_odd)

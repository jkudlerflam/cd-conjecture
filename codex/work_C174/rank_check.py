"""C174 check: period-p cells of a nearest-neighbour one-qubit-per-site chain.

Only one bond crosses each cell boundary, so the inter-cell Bloch block is
supported on 2 rows x 2 columns (two Majoranas, or a and a^dagger, per site).
Claim: det h(K) is a trigonometric polynomial of degree <= 2 in K for every
period p, fermionic (parity-class HFB) or bosonic (linear spin waves about a
period-p product reference). Hence at most 4 zeros counted with multiplicity:
c_MF <= 1 for fermions, and at most two soft momenta (+-k0, or {0, pi}) for
bosons, for any period. Checked here by the Fourier content of det h(K).
"""
import numpy as np
rng = np.random.default_rng(174)

def fermion_det(p, K):
    n = 2 * p
    A0 = np.zeros((n, n))
    # on-site and intra-cell nearest-neighbour Majorana bilinears only
    for s in range(p):
        for a in range(2):
            for b in range(2):
                i, j = 2 * s + a, 2 * s + b
                if i < j:
                    A0[i, j] = rng.normal()
                if s + 1 < p:
                    A0[2 * s + a, 2 * (s + 1) + b] = rng.normal()
    A0 = A0 - A0.T
    A1 = np.zeros((n, n))
    A1[n - 2:, :2] = rng.normal(size=(2, 2))   # site p of cell m to site 1 of cell m+1
    out = []
    for k in K:
        M = A0 + A1 * np.exp(1j * k) - A1.T * np.exp(-1j * k)
        out.append(np.linalg.det(1j * M))
    return np.array(out)

def boson_det(p, K):
    # bosonic BdG in basis (a_1..a_p, a_1^dag..a_p^dag); intra-cell NN + on-site,
    # inter-cell only between site p of cell m and site 1 of cell m+1
    n = 2 * p
    H0 = np.zeros((n, n), complex)
    for s in range(p):
        H0[s, s] = 5 + rng.random()           # keep it stable
        H0[p + s, p + s] = H0[s, s]
        g = rng.normal() + 1j * rng.normal()
        H0[s, p + s] = g; H0[p + s, s] = np.conj(g)
        if s + 1 < p:
            t = rng.normal() + 1j * rng.normal(); d = rng.normal() + 1j * rng.normal()
            H0[s, s + 1] += t; H0[s + 1, s] += np.conj(t)
            H0[p + s, p + s + 1] += np.conj(t); H0[p + s + 1, p + s] += t
            H0[s, p + s + 1] += d; H0[p + s + 1, s] += np.conj(d)
            H0[s + 1, p + s] += d; H0[p + s, s + 1] += np.conj(d)
    T = np.zeros((n, n), complex)
    t = rng.normal() + 1j * rng.normal(); d = rng.normal() + 1j * rng.normal()
    T[p - 1, 0] = t; T[2 * p - 1, p] = np.conj(t)
    T[p - 1, p] = d; T[2 * p - 1, 0] = np.conj(d)
    out = []
    for k in K:
        H = H0 + T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k)
        out.append(np.linalg.det(H))
    return np.array(out)

N = 64
K = 2 * np.pi * np.arange(N) / N
for name, f in (("fermion", fermion_det), ("boson", boson_det)):
    for p in (1, 2, 3, 4, 5, 7):
        worst = 0.0
        for trial in range(200):
            d = f(p, K)
            c = np.fft.fft(d) / N
            low = np.abs(c[[0, 1, 2, N - 1, N - 2]]).max()
            high = np.abs(np.delete(c, [0, 1, 2, N - 1, N - 2])).max()
            worst = max(worst, high / max(low, 1e-300))
        print(f"{name:8s} p={p}: max |harmonic n>2| / max |harmonic n<=2| over 200 models = {worst:.1e}")

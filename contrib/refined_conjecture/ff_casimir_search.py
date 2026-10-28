"""Free-fermion test of the Casimir bound sum_a c_a v_a <= C ||h||.
Chain of n orbitals per site (d = 2^n), Bloch Hamiltonian H(k) = A + B e^{ik} + B^dag e^{-ik},
bond h = c1^dag B c2 + h.c. + (1/2)(c1^dag A c1 + c2^dag A c2) (chemical potential inside A).
Many-body half spread of h = (1/2) * trace norm of the 2n x 2n single-particle bond matrix.
Each Fermi point contributes c = 1/2 with its velocity, so sum c_a v_a = (1/2) sum_FP |v_F|.
Random search plus local optimization of the ratio (sum c v)/||h||.
"""
import numpy as np
from scipy.optimize import minimize

K = 2048
ks = (np.arange(K) + 0.5) * 2 * np.pi / K


def unpack(x, n):
    a = x[: n * n].reshape(n, n); b = x[n * n: 2 * n * n].reshape(n, n) + 1j * x[2 * n * n:].reshape(n, n)
    A = (a + a.T) / 2
    return A, b


def ratio(x, n):
    A, B = unpack(x, n)
    Hk = A[None] + B[None] * np.exp(1j * ks)[:, None, None] + B.conj().T[None] * np.exp(-1j * ks)[:, None, None]
    e = np.linalg.eigvalsh(Hk)  # (K, n), sorted
    vsum = 0.0
    for b in range(n):
        eb = e[:, b]
        s = np.sign(eb); idx = np.nonzero(s != np.roll(s, -1))[0]
        for i in idx:
            j = (i + 1) % K
            vsum += abs(eb[j] - eb[i]) / (2 * np.pi / K)
    M = np.block([[A / 2, B], [B.conj().T, A / 2]])
    hnorm = 0.5 * np.abs(np.linalg.eigvalsh(M)).sum()
    return 0.5 * vsum / hnorm


import sys
ns = [int(a) for a in sys.argv[1:]] or [1, 2, 3, 4]
rng = np.random.default_rng(1)
for n in ns:
    best = 0
    for trial in range(15):
        x0 = rng.normal(size=3 * n * n)
        r0 = ratio(x0, n)
        res = minimize(lambda x: -ratio(x, n), x0, method="Nelder-Mead", options=dict(maxiter=800, xatol=1e-6, fatol=1e-8))
        best = max(best, r0, -res.fun)
    print(f"n={n} (d={2**n}): max found (sum c v)/||h|| = {best:.4f}   [2 pi = {2*np.pi:.4f}]", flush=True)

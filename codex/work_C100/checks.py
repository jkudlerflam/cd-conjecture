"""Small finite-volume and quadrature checks of the C100 derivations.

These are normalization checks, not evidence for an infrared identification.
Run with python3 checks.py. Only numpy and scipy are used for computation.
"""
from pathlib import Path
import numpy as np
from scipy.integrate import quad


def comm(a, b):
    return a @ b - b @ a


def embed_bond(h, i, n):
    out = np.zeros((2**n, 2**n), complex)
    j = (i + 1) % n
    for col in range(2**n):
        bits = [(col >> k) & 1 for k in range(n)]
        local = 2 * bits[i] + bits[j]
        for a in range(2):
            for b in range(2):
                row = col + ((a - bits[i]) << i) + ((b - bits[j]) << j)
                out[row, col] += h[2 * a + b, local]
    return out


rng = np.random.default_rng(100)
n = 6
max_cont = max_spec = max_bound = max_mid = 0.0
min_moment = np.inf
for trial in range(4):
    raw = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    h = (raw + raw.conj().T) / 2
    ev = np.linalg.eigvalsh(h)
    h = (h - (ev[-1] + ev[0]) * np.eye(4) / 2) / ((ev[-1] - ev[0]) / 2)
    hs = [embed_bond(h, i, n) for i in range(n)]
    H = sum(hs)
    js = [1j * comm(hs[(i - 1) % n], hs[i]) for i in range(n)]
    for i in range(n):
        err = np.linalg.norm(1j * comm(H, hs[i]) - js[i] + js[(i + 1) % n])
        max_cont = max(max_cont, err)
    e, U = np.linalg.eigh(H)
    omega = e[:, None] - e[None, :]
    for beta in (0.0, 0.3, 2.0, 8.0):
        p = np.exp(-beta * (e - e[0]))
        p /= p.sum()
        rho = (U * p) @ U.conj().T
        for k in (1, 2, 3):
            q = 2 * np.pi * k / n
            hq = sum(np.exp(-1j * q * i) * hs[i] for i in range(n))
            direct = np.trace(rho @ comm(hq.conj().T, comm(H, hq))).real / (2 * n)
            aq, am = U.conj().T @ hq @ U, U.conj().T @ hq.conj().T @ U
            spectral = np.sum(omega * (np.abs(aq)**2 + np.abs(am)**2) * p[None, :]).real / (2 * n)
            max_spec = max(max_spec, abs(direct - spectral))
            max_bound = max(max_bound, direct / (8 * q**2))
            min_moment = min(min_moment, direct)
        jj = U.conj().T @ js[0] @ U
        midpoint = np.sum(np.sqrt(p[:, None] * p[None, :]) * np.abs(jj)**2).real
        max_mid = max(max_mid, midpoint / 4)

# The finite-temperature chiral spectrum follows from the midpoint sech^4 kernel.
c, v, temp = 1.7, 2.3, 0.41
midspace = quad(lambda x: c * np.pi**2 * temp**4 / 4
                / np.cosh(np.pi * temp * x / v)**4, -60, 60)[0]
mid_expected = np.pi * c * v * temp**3 / 3
qs = [0.01, 0.4, 1.2]
fourier_errors = []
for q in qs:
    val = quad(lambda x: np.cos(q * x) * c * np.pi**2 * temp**4 / 4
               / np.cosh(np.pi * temp * x / v)**4, -60, 60)[0]
    K = c * v**2 / (24 * np.pi) * q * (q*q + (2*np.pi*temp/v)**2)
    expected = v**2 * K / np.sinh(v*q/(2*temp))
    fourier_errors.append(abs(val - expected))

assert max_cont < 1e-11
assert max_spec < 1e-11
assert min_moment > -1e-11
assert max_bound <= 1 + 1e-12
assert max_mid <= 1 + 1e-12
assert abs(midspace - mid_expected) < 1e-11
assert max(fourier_errors) < 1e-11

text = f'''# C100 checks

[Numerical checks, not proofs of CFT matching.] Four seeded random Hermitian
two-qubit bonds, each centered and normalized to operator norm one, on six-site
rings. Temperatures sampled via beta = 0, 0.3, 2, 8. Three nonzero momenta.
The centered norm is used, not a numerical gauge optimization for M.

- Maximum continuity-equation matrix residual: {max_cont:.3e}.
- Maximum first-moment discrepancy, double commutator versus spectral sum: {max_spec:.3e}.
- Smallest first moment: {min_moment:.3e} (roundoff allowed at beta = 0).
- Largest ratio m1/(8 m^3 q^2): {max_bound:.6f}, with m = 1.
- Largest local midpoint correlator divided by 4 m^4: {max_mid:.6f}.
- CFT midpoint spatial integral error: {abs(midspace-mid_expected):.3e}.
- Largest CFT midpoint Fourier-transform error: {max(fourier_errors):.3e}.

All assertions passed. No large computation, infrared fit, or Drude-limit
extrapolation was performed. Analytic checks of scaling, chiral conventions,
and the absence of a temperature-independent conclusion are in PROOF.md.
'''
Path('CHECKS.md').write_text(text)
print(text)

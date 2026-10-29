"""Small analytic checks for C107. No many-body computation above 256 states."""
import numpy as np
from scipy.integrate import quad
from pathlib import Path


def xy_energy(L, gamma):
    k = 2*np.pi*(np.arange(L)+0.5)/L
    return -np.sum(np.sqrt((1-np.cos(k))**2 + gamma**2*np.sin(k)**2))


def xy_bulk(gamma):
    return -quad(lambda k: np.sqrt((1-np.cos(k))**2 + gamma**2*np.sin(k)**2),
                 0, np.pi, epsabs=2e-13, epsrel=2e-13)[0]/np.pi


def spin_energy(L, gamma):
    # H = -sum [(1+gamma)XX/2 + (1-gamma)YY/2 + Z].
    n = 1 << L
    H = np.zeros((n,n))
    for s in range(n):
        H[s,s] = -(L-2*s.bit_count())
        for j in range(L):
            k = (j+1) % L
            a, b = (s >> j) & 1, (s >> k) & 1
            H[s ^ (1 << j) ^ (1 << k),s] -= 1.0 if a != b else gamma
    return np.linalg.eigvalsh(H)[0]


out = ['# C107 checks', '', '[Numerical] Small spectral and quadrature checks of the analytic derivations. Run with `python3 checks.py`. These are not extrapolation evidence for the general conjecture.', '']
errors = []
for gamma in [0.25, 0.5, 1.0]:
    for L in [4,6,8]:
        errors.append(abs(spin_energy(L,gamma)-xy_energy(L,gamma)))
assert max(errors) < 1e-11
out += [f'[Numerical] Nine full spin matrices, L = 4, 6, 8 and gamma = 0.25, 0.5, 1, agree with the antiperiodic fermion sum. Largest dimension 256, largest energy discrepancy {max(errors):.3g}.', '', '| gamma | L | L^2 (Q_L - gamma) | predicted coefficient |', '|:--|--:|--:|--:|']
for gamma in [0.1, 0.25, 0.5, 1.0]:
    e = xy_bulk(gamma)
    predicted = 7*np.pi**2*(4*gamma**2-3)/(240*gamma)
    qs = []
    for L in [64,128,256]:
        Q = 6/np.pi*L*(L*e-xy_energy(L,gamma))
        qs.append(Q)
        out.append(f'| {gamma:g} | {L} | {L**2*(Q-gamma):.9f} | {predicted:.9f} |')
    assert (qs[-1]-qs[-2])*predicted < 0
    assert abs((256**2*(qs[-1]-gamma)-predicted)/predicted) < 0.03
out += ['', '[Numerical] The sign and large-L coefficient agree with the Euler-Maclaurin derivation for all four anisotropies. The smallest anisotropy has the longest crossover.', '', '| L | half-filled ferromagnetic XX Q_L, t = 1 |', '|--:|--:|']
for L in [63,64,127,128]:
    E = -2/np.sin(np.pi/L) if L % 2 == 0 else -1/np.sin(np.pi/(2*L))
    Q = 6/np.pi*L*(-2*L/np.pi-E)
    out.append(f'| {L} | {Q:.9f} |')
out += ['', '[Numerical] Even and odd XX lengths approach 2 and 1/2 respectively, in agreement with the exact formulas.', '']
for L in [1,2,10,100,1000]:
    t = np.linspace(0,1,1001)
    e2 = -6*L+2*L*t*(3-L*t*t)/(1+L*t*t)**3
    slope = -6*L*t-(1-L*t*t)/(1+L*t*t)**2
    assert e2.max() <= 1e-12
    assert np.max(np.abs(slope)) <= 7*L+1e-12
for n in [1,2,10,100]:
    for m in [1,3,17,100]:
        t = np.linspace(0,1,1001)
        d = lambda L: t/(1+L*t*t)
        defect = d(n)+d(m)-d(n+m)
        assert defect.min() >= -1e-12
        assert defect.max() <= 2+1e-12
out += ['[Numerical] The abstract countermodel derivative and additivity checks pass on the stated grids. Their proofs do not rely on sampling.', '', '[Proved scope] No large diagonalization, finite-size fit, or cluster job was used. No physical bond family with unbounded A was constructed.', '']
Path('CHECKS.md').write_text('\n'.join(out))
print('All checks passed. Maximum spin/fermion discrepancy:', max(errors))

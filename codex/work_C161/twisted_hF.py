"""C161 (Claude Opus agent, 2026-11-08). The singular set through h_F is a curve: J = -R(phi) (R a rotation by phi about z)
is h_F after the site-dependent rotation sigma_j -> R(j phi) sigma_j. Check: open-chain spectra of J = -R(phi), B = 0,
against h_F, L = 8. Also: XXZ Luttinger parameter thresholds of C159 Proposition C in terms of Delta_AF and eps.
"""
import numpy as np
from functools import reduce
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0, -1.0]).astype(complex)
P = [sx, sy, sz]
def op(o, j, L):
    return reduce(np.kron, [o if i == j else np.eye(2) for i in range(L)])
def H(J, L):
    S = [[op(P[a], j, L) for a in range(3)] for j in range(L)]
    return sum(J[a, b] * S[j][a] @ S[j + 1][b] for j in range(L - 1) for a in range(3) for b in range(3))
L = 8
e0 = np.linalg.eigvalsh(H(-np.eye(3), L))
for phi in (0.3, 1.0, 2.0, np.pi):
    c, s = np.cos(phi), np.sin(phi)
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    e = np.linalg.eigvalsh(H(-R, L))
    print(f'phi={phi:.3f}  max|spec(-R(phi)) - spec(h_F)| = {np.abs(e - e0).max():.2e}   (J = -R(phi) has J_perp = {-c:+.3f}, D_z = {s:+.3f}, J_z = -1)')
# C159 thresholds: K = pi / (2 (pi - arccos Delta_AF)) for the XXZ chain (K = 1 at Delta = 0, K -> inf at Delta -> -1)
for K in (1.0, 15/8, 3.0, 5.0):
    g = np.pi - np.pi / (2 * K); D = np.cos(g)
    eps_dm = np.sqrt(1 / D ** 2 - 1)        # pure-DM line of C160: Delta_AF = -1/sqrt(1 + eps^2)
    print(f'K={K:.3f}: Delta_AF={D:+.4f}, eta=1+Delta={1+D:.4f}, pure-DM eps={eps_dm:.3f}, hbar_eff(pure DM)=2 eps={2*eps_dm:.3f}, K*hbar_eff={K*2*eps_dm:.2f}')

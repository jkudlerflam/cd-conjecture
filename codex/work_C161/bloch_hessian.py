"""C161 (Claude Opus agent, 2026-11-08). Test of C160 Lemma M2 on a non-planar texture with DM along the moment:
the tilted (conical-like) CSL of the K-stratum functional W = -h n_x + (kappa/2) n_z^2 + (mu/2) n_y^2, D = z.
Bloch Hessian bands lambda_0(k), lambda_1(k) of the period-optimized texture (tangent two-component, gyroscopic DM).
M2 predicts: lambda_0(0) = 0 (phason), lambda_0(k) > 0 for 0 < k <= pi/T, minimum at k = 0.
"""
import sys
import numpy as np
from scipy.optimize import minimize_scalar
from soliton_pair import ring_min, ring_energy

h, kappa, mu = (float(a) for a in sys.argv[1:4]) if len(sys.argv) > 3 else (0.015, 0.05, 2.3)
dx0 = 0.1
res = minimize_scalar(lambda T: ring_min(h, kappa, mu, T, 1, [0.6], dx=dx0)[0] / T, bounds=(7.0, 12.0),
                      method='bounded', options=dict(xatol=1e-3))
T = res.x
E, n = ring_min(h, kappa, mu, T, 1, [0.6], dx=dx0)
N = len(n); dx = T / N
print(f'h={h} kappa={kappa} mu={mu}: T={T:.4f}, e={E/T:+.6f}, max|n_z|={np.abs(n[:,2]).max():.3f}, <n_z>={n[:,2].mean():+.3f}, max|D.n|={np.abs(n[:,2]).max():.3f}')
a = np.cross(n, [0.0, 0.0, 1.0]); a[np.linalg.norm(a, axis=1) < 1e-8] = [1, 0, 0]
e1 = a / np.linalg.norm(a, axis=1)[:, None]; e2 = np.cross(n, e1)

def g(ab):
    A = ab.reshape(N, 2)
    v = n + A[:, :1] * e1 + A[:, 1:] * e2
    G = ring_energy(v.ravel(), dx, h, kappa, mu)[1].reshape(N, 3)
    return np.stack([np.sum(G * e1, 1), np.sum(G * e2, 1)], 1).ravel()

d = 1e-5
H = np.zeros((2 * N, 2 * N))
for j in range(2 * N):
    ej = np.zeros(2 * N); ej[j] = d
    H[:, j] = (g(ej) - g(-ej)) / (2 * d)
H = 0.5 * (H + H.T)
print('residual gradient norm', np.linalg.norm(g(np.zeros(2 * N))))
# Bloch: couplings between the last point block and the first point block get phase e^{ikT}
B = np.zeros((2 * N, 2 * N), bool)
B[2 * (N - 1):, :2] = True
for k in np.linspace(0, np.pi / T, 9):
    Hk = H.astype(complex)
    Hk[2 * (N - 1):, :2] *= np.exp(1j * k * T)
    Hk[:2, 2 * (N - 1):] *= np.exp(-1j * k * T)
    w = np.linalg.eigvalsh(Hk) / dx
    print(f'k*T/pi={k*T/np.pi:.3f}  lambda_0={w[0]:+.6f}  lambda_1={w[1]:+.6f}')

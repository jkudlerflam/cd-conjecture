"""Random-instance tests of two finite inequalities used by C44/C45.

1. C44 PROOFS (13): Schur-complement equalization bound
   ||s Z^dag (A+sW)^-1 Z - z^dag B0^-1 z|| <= eps_eq, with A = (H_ch - e) x I
   (charge ground chi, gap delta), W = sum_j D_j x h_j >= 0 (D_j 0/1 adjacency,
   h_j >= 0), B = PWP = sum_j d_j h_j, B0 = sum_j h_j, Z = chi x z + Y.
   Built in the same block structure as the wrapper, with random data.
2. C45 PROOF (22): if ||K0^-1/2 (K-K0) K0^-1/2|| <= 1/4 then the Gaussian
   dimension metric G(K) = K^-1 # (S^T K S) obeys G(K) >= sqrt(3/5) G(K0),
   for stiffness matrices with velocity ratios up to 1e4.
"""
import json
import numpy as np
from scipy.linalg import sqrtm

rng = np.random.default_rng(2026)


def psd(n, scale=1.0, rank=None):
    r = rank or n
    m = rng.normal(size=(n, r))
    return scale * m @ m.T / r


def equalize_trial(mch=6, c=5, nbond=3, nz=2):
    Hc = rng.normal(size=(mch, mch)); Hc = (Hc + Hc.T) / 2
    ev, evec = np.linalg.eigh(Hc)
    ev = ev - ev[0]
    ev[1:] = np.maximum(ev[1:], rng.uniform(0.05, 1.0))   # gap delta
    chi = evec[:, 0]
    Hc = evec @ np.diag(ev) @ evec.T
    delta = ev[1]
    A = np.kron(Hc, np.eye(c))
    Ds = [np.diag((rng.uniform(size=mch) > rng.uniform(0, 0.3)).astype(float)) for _ in range(nbond)]
    hs = [psd(c, rng.uniform(0.2, 2.0), rank=rng.integers(1, c + 1)) for _ in range(nbond)]
    hs[0] = hs[0] + np.eye(c) * rng.uniform(0.3, 1.0)       # pair core gap
    W = sum(np.kron(D, h) for D, h in zip(Ds, hs))
    d = [chi @ D @ chi for D in Ds]
    B = sum(dj * h for dj, h in zip(d, hs))
    B0 = sum(hs)
    a = min(np.linalg.eigvalsh(B)[0], np.linalg.eigvalsh(B0)[0])
    zeta2 = max(1 - dj for dj in d)
    Lb = sum(np.linalg.norm(h, 2) for h in hs)
    Lc = np.linalg.norm(W, 2)
    z = rng.normal(size=(c, nz)) * rng.uniform(0.1, 2)
    Pch = np.outer(chi, chi)
    Yc = rng.normal(size=(mch * c, nz)) * rng.uniform(0, 3)
    Y = np.kron(np.eye(mch) - Pch, np.eye(c)) @ Yc
    Z = np.kron(chi[:, None], z) + Y
    x, y = np.linalg.norm(z, 2), np.linalg.norm(Y, 2)
    out = []
    for frac in (1e-3, 0.05, 0.3, 1.0):
        t = frac * a / (2 * Lc ** 2)
        s = t * delta
        lhs = np.linalg.norm(s * Z.T @ np.linalg.solve(A + s * W, Z) - z.T @ np.linalg.solve(B0, z), 2)
        rhs = (x ** 2 * zeta2 * Lb / a ** 2 + 2 * t * Lc ** 2 * x ** 2 / a ** 2
               + (2 / a) * (2 * t * Lc * x * y + t ** 2 * Lc ** 2 * y ** 2) + t * y ** 2)
        out.append(lhs / rhs)
    return out


ratios = np.array([equalize_trial() for _ in range(3000)])
res = dict(equalize_trials=int(ratios.size), equalize_max_ratio=float(ratios.max()),
           equalize_median_ratio=float(np.median(ratios)))


def gmean(X, Y):
    Xh = np.real(sqrtm(X)); Xmh = np.linalg.inv(Xh)
    return Xh @ np.real(sqrtm(Xmh @ Y @ Xmh)) @ Xh


worst = 1e9
for _ in range(2000):
    m = rng.integers(1, 4)
    v = 10 ** rng.uniform(-4, 0, size=m)
    Kl = rng.uniform(0.5, 2, size=m)
    O = np.linalg.qr(rng.normal(size=(m, m)))[0]
    A0 = O @ np.diag(v / Kl) @ O.T; B0 = O @ np.diag(v * Kl) @ O.T
    K0 = np.block([[A0, 0.1 * np.sqrt(np.outer(v, v)) * rng.normal(size=(m, m))], [np.zeros((m, m)), B0]])
    K0 = (K0 + K0.T) / 2
    if np.linalg.eigvalsh(K0)[0] <= 0:
        continue
    E = rng.normal(size=(2 * m, 2 * m)); E = (E + E.T) / 2; E *= 0.25 / np.linalg.norm(E, 2)
    K0h = np.real(sqrtm(K0))
    K = K0h @ (np.eye(2 * m) + E) @ K0h
    S = np.block([[np.zeros((m, m)), np.eye(m)], [-np.eye(m), np.zeros((m, m))]])
    G0 = gmean(np.linalg.inv(K0), S.T @ K0 @ S)
    G = gmean(np.linalg.inv(K), S.T @ K @ S)
    G0mh = np.linalg.inv(np.real(sqrtm(G0)))
    worst = min(worst, np.linalg.eigvalsh(G0mh @ G @ G0mh)[0])
res.update(gauss_min_ratio=float(worst), sqrt_3_5=float(np.sqrt(0.6)),
           min_dimension_of_D3_vertex=float(3 * worst))
json.dump(res, open('ref_equalize_gauss.json', 'w'), indent=1)
print(json.dumps(res, indent=1))

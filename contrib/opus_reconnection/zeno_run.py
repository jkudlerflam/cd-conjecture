"""Continuous-pinching (Zeno) limit of G at fixed beta, any chain.

Proved in README Section 4: as N -> infinity at fixed beta,
  S_pin -> sum_s z_ss^2,  Q_D -> sum_{s,t} z_st z_ts,
with z_st = Tr exp(-beta H_open(s,t)) the open (L-1)-site chain with the
boundary fields produced by pinned seam values s (left) and t (right).
Hence G_inf <= 0. This script checks the convergence of G_N to G_inf, the
sign of G_inf, and the growth of the pinching loss eps_L ~ beta (E_Zeno-E_L).
Output: zeno_results.json."""
import json
from itertools import product
from pathlib import Path
import numpy as np
from scipy.linalg import expm, eigh

ROOT = Path(__file__).resolve().parent

def bond(model, delta=1.0):
    if model == 'su3':
        P = np.zeros((9, 9))
        for a in range(3):
            for b in range(3): P[3 * b + a, 3 * a + b] = 1
        return P, 3
    sx = np.array([[0, 1], [1, 0]]) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2
    sz = np.diag([.5, -.5])
    h = np.kron(sx, sx) + np.kron(sy, sy) + delta * np.kron(sz, sz)
    return h.real, 2

def ring(h, d, n):
    H = np.zeros((d ** n, d ** n))
    for i in range(n):
        j = (i + 1) % n
        # operator on sites i,j (site 0 is the most significant factor)
        perm = list(range(n)); 
        op = h.reshape(d, d, d, d)
        T = np.zeros((d,) * (2 * n))
        H += _embed(op, i, j, n, d)
    return H

def _embed(op, i, j, n, d):
    dim = d ** n
    idx = np.indices((d,) * n).reshape(n, -1).T
    M = np.zeros((dim, dim))
    for col, conf in enumerate(idx):
        a, b = conf[i], conf[j]
        for a2 in range(d):
            for b2 in range(d):
                v = op[a2, b2, a, b]
                if v != 0:
                    c2 = conf.copy(); c2[i] = a2; c2[j] = b2
                    row = np.ravel_multi_index(c2, (d,) * n)
                    M[row, col] += v
    return M

def seam_labels(n, d, sites):
    idx = np.indices((d,) * n).reshape(n, -1).T
    return sum(d ** k * idx[:, s] for k, s in enumerate(sites))

def G_finite(HL, HD, d, L, beta, N):
    lab = seam_labels(L, d, [0]); labD = seam_labels(2 * L, d, [0, L])
    tau = beta / N
    shift = lambda H: H - eigh(H, eigvals_only=True)[0] * np.eye(len(H))
    F = expm(-tau * shift(HL) / 4); E = F @ F
    q = []
    for hist in product(range(d), repeat=N):
        w = F.copy()
        for j, s in enumerate(hist):
            w = (lab == s)[:, None] * w
            w = (F if j == N - 1 else E) @ w
        q.append(np.sum(w * w))
    q = np.array(q)
    eD = expm(-tau * shift(HD) / 2); mask = labD[:, None] == labD[None, :]
    rho = eD.copy()
    for j in range(N):
        rho *= mask
        if j < N - 1: rho = eD @ rho @ eD
    QD = np.sum(rho * eD.T)
    EL = eigh(HL, eigvals_only=True)[0]; ED_ = eigh(HD, eigvals_only=True)[0]
    zL = np.sum(np.exp(-beta * (eigh(HL, eigvals_only=True) - EL)))
    return float(beta * (ED_ - 2 * EL) + np.log(q @ q) - np.log(QD)), float(np.log(zL) - np.log(q.sum()))

def G_zeno(HL, HD, d, L, beta):
    lab = seam_labels(L, d, [0]); labD = seam_labels(2 * L, d, [0, L])
    EL = eigh(HL, eigvals_only=True)[0]; ED_ = eigh(HD, eigvals_only=True)[0]
    zss = []
    for s in range(d):
        k = lab == s
        zss.append(np.sum(np.exp(-beta * (eigh(HL[np.ix_(k, k)], eigvals_only=True) - EL))))
    zD = 0.0
    for st in range(d * d):
        k = labD == st
        zD += np.sum(np.exp(-beta * (eigh(HD[np.ix_(k, k)], eigvals_only=True) - ED_)))
    zss = np.array(zss)
    # Zeno ground energy shift of the pinched L ring
    Ez = min(eigh(HL[np.ix_(lab == s, lab == s)], eigvals_only=True)[0] for s in range(d))
    return float(beta * (ED_ - 2 * EL) + np.log(zss @ zss) - np.log(zD)), float(Ez - EL)

def main():
    rows = []
    for name, model, delta, L, M in [('XX', 'xxz', 0.0, 4, .5), ('Heisenberg', 'xxz', 1.0, 4, .5),
                                     ('SU3', 'su3', 1.0, 3, 1.)]:
        h, d = bond(model, delta)
        HL, HD = ring(h, d, L), ring(h, d, 2 * L)
        for beta in [2 * L / M * 0.25, L / M]:   # tau-independent fixed beta
            Ginf, dEz = G_zeno(HL, HD, d, L, beta)
            Ns = range(1, 11) if d == 2 else range(1, 8)
            seq = [G_finite(HL, HD, d, L, beta, N) for N in Ns]
            rows.append(dict(model=name, L=L, beta=beta, G_inf=Ginf, Ezeno_minus_EL=dEz,
                             N=list(Ns), G_N=[s[0] for s in seq], epsL_N=[s[1] for s in seq]))
            print(name, 'beta', beta, 'G_inf', round(Ginf, 5), 'E_Z-E_L', round(dEz, 5))
            print('   G_N  ', [round(s[0], 4) for s in seq])
            print('   epsL ', [round(s[1], 4) for s in seq])
            assert Ginf <= 1e-12
    (ROOT / 'zeno_results.json').write_text(json.dumps(rows, indent=1) + '\n')

if __name__ == '__main__':
    main()

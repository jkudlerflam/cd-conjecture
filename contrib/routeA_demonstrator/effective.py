#!/usr/bin/env python3
"""Second-order effective logical bond of two demonstrator cells, and the exact two-cell check.

H(lam) = H0 (x) 1 + 1 (x) H0 + lam V,  V on (pR of cell 1, pL of cell 2).
Second order: H2 = - P V Q (H0 (x) 1 + 1 (x) H0)^{-1} Q V P on the 4-dim code (E0 = 0).
Decomposition on the logical pair:  H2 = a0 + (Jxy/2)(XX + YY) + Jzz ZZ + h (Z1 + 1Z) + g (Z1 - 1Z).
Exact check: des Cloizeaux effective Hamiltonian from the 4 lowest exact two-cell eigenvectors.
Claude (Opus 5.5), 2026-10-27.
"""
import sys, json, argparse
import numpy as np
from scipy import linalg
import cell as C

PAU = dict(I=np.eye(2), X=C.X, Y=np.array([[0, -1j], [1j, 0]]), Z=C.Z)


def decompose(M):
    """Pauli coefficients of a 4x4 logical operator (cell 1 first)."""
    out = {}
    for a in 'IXYZ':
        for b in 'IXYZ':
            out[a + b] = float(np.real(np.trace(np.kron(PAU[a], PAU[b]).conj().T @ M)) / 4)
    return out


def second_order(D0=1.0, D1=1.0, kap=0.5, Jp=1.0, Jz=1.0):
    H00, V, code, qc = C.inter_cell(lam=1.0, Jp=Jp, Jz=Jz, D0=D0, D1=D1, kap=kap)
    G = np.kron(code, code)                          # 1024 x 4
    Pc = G @ G.T
    src = V @ G                                      # first order would be G^T src
    first = G.T @ src
    src = src - Pc @ src
    # solve H00 x = src on the complement (H00 >= gap on Q, and src is in Q)
    w, U = linalg.eigh(H00)
    mask = w > 1e-10
    Ui = U[:, mask]
    x = Ui @ ((Ui.T @ src) / w[mask, None])
    H2 = -src.T @ x
    return first, H2, w[mask].min()


def exact_eff(lam, D0=1.0, D1=1.0, kap=0.5, Jp=1.0, Jz=1.0):
    H00, V, code, qc = C.inter_cell(lam=1.0, Jp=Jp, Jz=Jz, D0=D0, D1=D1, kap=kap)
    G = np.kron(code, code)
    H = H00 + lam * V
    w, U = linalg.eigh(H)
    low = U[:, :4]
    Phi = G.T @ low                                  # 4x4 overlaps
    s = linalg.svdvals(Phi)
    # Loewdin / des Cloizeaux:  O = Phi (Phi^T Phi)^{-1/2}, Heff = O diag(E) O^T
    Sm = linalg.inv(linalg.sqrtm(Phi.T @ Phi))
    O = Phi @ Sm
    Heff = O @ np.diag(w[:4]) @ O.T
    return np.real(Heff), w[:6], s.min()


def summary(H2):
    d = decompose(H2)
    return dict(a0=d['II'], Jxy=d['XX'] * 2, Jxy_check=d['YY'] * 2, Jzz=d['ZZ'],
                h=(d['ZI'] + d['IZ']) / 2, g=(d['ZI'] - d['IZ']) / 2,
                other=float(np.sqrt(sum(v ** 2 for k, v in d.items() if k not in ('II', 'XX', 'YY', 'ZZ', 'ZI', 'IZ')))),
                anis=d['XX'] - d['YY'])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--D0', type=float, default=1.0)
    ap.add_argument('--D1', type=float, default=1.0)
    ap.add_argument('--kap', type=float, default=0.5)
    ap.add_argument('--Jp', type=float, default=1.0)
    ap.add_argument('--Jz', type=float, default=1.0)
    ap.add_argument('--lams', type=float, nargs='*', default=[0.4, 0.2, 0.1, 0.05, 0.025])
    a = ap.parse_args()
    first, H2, gap = second_order(a.D0, a.D1, a.kap, a.Jp, a.Jz)
    print('cell-pair gap', gap, ' |first order| =', np.abs(first).max())
    s2 = summary(H2)
    print('second order:', {k: round(v, 8) for k, v in s2.items()})
    print('Delta_xxz = Jzz/(Jxy/2) =', s2['Jzz'] / (s2['Jxy'] / 2))
    for lam in a.lams:
        Heff, w, smin = exact_eff(lam, a.D0, a.D1, a.kap, a.Jp, a.Jz)
        r = Heff - lam ** 2 * H2
        print(f'lam={lam:6.3f}  low6={np.round(w, 6)}  min overlap sv={smin:.4f}  '
              f'|Heff - lam^2 H2|={np.linalg.norm(r, 2):.3e}  /lam^3={np.linalg.norm(r, 2) / lam ** 3:.4f}  '
              f'/lam^4={np.linalg.norm(r, 2) / lam ** 4:.4f}')

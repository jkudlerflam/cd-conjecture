#!/usr/bin/env python3
"""Cell selection by energetics: run-unit energies of the demonstrator chain.

Site types are conserved exactly (no term changes P <-> B), so a ring is a sequence of B's separated
by runs of k >= 0 P's (k = 0 means a BB bond).  Writing H as a sum over run units
  R_k = H_R(sR, c ; p_1) + lam sum_i V(p_i, p_i+1) + H_L(p_k ; sL', c')      (k >= 1)
  R_0 = K  (a BB bond; its sR, sL are unchecked)
plus -mu per B gives the exact operator inequality  H >= sum_B (-mu) + sum_runs e_k,  e_k = min spec R_k,
because R_k >= e_k as operators even though neighbouring units share the logical qubit c of a B.
On a ring every site belongs to one B or one run, so the energy per site of any configuration is at
least a weighted average of u_k = (-mu + e_k)/(k + 1).  The period-3 pattern (all k = 2) is selected
when u_2 < u_k for every k != 2, with margins printed below; the true period-3 energy lies above
u_2 by the logical-chain energy, which is O(lam^2) and is printed for comparison (from a 2-cell ring
bound).  Claude (Opus 5.5), 2026-10-27.
"""
import argparse, json
import numpy as np
from scipy.sparse.linalg import eigsh
import cell as C


def run_operator(k, lam, D0, D1, kap, Jp, Jz):
    """R_k on (sR, c, p_1..p_k, sL', c'), qubit order exactly that, dimension 2^(k+4)."""
    n = k + 4
    hp = C.h_pair(D0, D1, kap)                       # (p, s, c)
    I = lambda m: np.eye(2 ** m)
    # right check of the left B: p = p_1 (pos 2), s = sR (0), c (1)
    HR = C.perm_op(np.kron(hp, I(n - 3)), n, [2, 0, 1] + list(range(3, n)))
    # left check of the right B: p = p_k (pos k+1), s = sL' (k+2), c' (k+3)
    rest = [i for i in range(n) if i not in (k + 1, k + 2, k + 3)]
    HL = C.perm_op(np.kron(hp, I(n - 3)), n, [k + 1, k + 2, k + 3] + rest)
    H = HR + HL
    Vb = Jp * (np.kron(C.SP, C.SM) + np.kron(C.SM, C.SP)) + Jz * np.kron(C.Z, C.Z)
    for i in range(k - 1):
        pos = 2 + i
        H += lam * C.perm_op(np.kron(Vb, I(n - 2)), n, [pos, pos + 1] + [j for j in range(n) if j not in (pos, pos + 1)])
    return H


def frustrated(D0, D1, kap):
    return np.linalg.eigvalsh(run_operator(1, 0.0, D0, D1, kap, 0, 0))[0]


def unit_energies(lam, D0, D1, kap, Jp, Jz, kmax=7):
    e = {}
    for k in range(1, kmax + 1):
        H = run_operator(k, lam, D0, D1, kap, Jp, Jz)
        e[k] = float(np.linalg.eigvalsh(H)[0]) if H.shape[0] <= 1024 else float(eigsh(H, k=1, which='SA')[0][0])
    return e


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--D0', type=float, default=1.0)
    ap.add_argument('--D1', type=float, default=1.0)
    ap.add_argument('--kap', type=float, default=0.6)
    ap.add_argument('--Jp', type=float, default=2.0)
    ap.add_argument('--Jz', type=float, default=-1.0)
    ap.add_argument('--mu', type=float, default=0.55)
    ap.add_argument('--K', type=float, default=2.0)
    ap.add_argument('--lams', type=float, nargs='+', default=[0.0, 0.05, 0.1, 0.15])
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    print('frustrated-port cost eps =', frustrated(a.D0, a.D1, a.kap))
    res = {}
    for lam in a.lams:
        e = unit_energies(lam, a.D0, a.D1, a.kap, a.Jp, a.Jz)
        u = {0: -a.mu + a.K}
        u.update({k: (-a.mu + e[k]) / (k + 1) for k in e})
        # tail k >= 8:  R_k >= [H_R + lam V(p1,p2)] + [lam V(p_k-1,p_k) + H_L] + lam sum_middle V >= 2f - lam (k-3)
        # since every middle V >= -lam (its spectrum is {-1,-1,-1,3}); f = min spec of one bracket.
        f = float(np.linalg.eigvalsh(run_operator(2, lam, a.D0, a.D1, a.kap, a.Jp, a.Jz)
                                     - C.perm_op(np.kron(C.h_pair(a.D0, a.D1, a.kap), np.eye(8)), 6, [3, 4, 5, 0, 1, 2]))[0])
        tail = min((-a.mu + 2 * f - lam * (k - 3)) / (k + 1) for k in range(8, 2000))
        u['tail>=8'] = tail
        margins = {k: u[k] - u[2] for k in u if k != 2}
        print(f'lam={lam}: e_k =', {k: round(v, 5) for k, v in e.items()})
        print('   u_k =', {k: round(v, 5) for k, v in u.items()}, ' min margin', round(min(margins.values()), 5),
              'at k =', min(margins, key=margins.get))
        res[lam] = dict(e=e, u={str(k): v for k, v in u.items()}, margins={str(k): v for k, v in margins.items()}, f=f)
    if a.out:
        json.dump(dict(args=vars(a), res={str(k): v for k, v in res.items()}), open(a.out, 'w'), indent=1)

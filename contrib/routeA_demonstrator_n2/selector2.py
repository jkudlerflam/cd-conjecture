#!/usr/bin/env python3
"""Cell selection for the 13-letter bond: run-unit energies, as in ../routeA_demonstrator/selector.py.

Run unit R_k (k >= 1 P letters between two B's) on (B_left logical, P_1..P_k, B_right logical):
  H_R(P_1 ; B_left) + lam sum_i V(x_i, x_i+1) + H_L(P_k ; B_right).
H >= sum_B (-mu) + sum_runs e_k exactly (e_k = min spec R_k), so the energy per site of any type
configuration is a weighted average of u_k = (-mu + e_k)/(k + 1) (u_0 = -mu + K for a BB bond).
Tail k >= 5: R_k >= 2 f - |v_min| lam (k - 3), with f the lowest eigenvalue of [H_R(P_1;B) + lam V(x_1,x_2)]
and v_min the lowest eigenvalue of V (= -2 for Jp = Jz = 1).
Claude (Opus 5.5), 2026-10-28.
"""
import argparse, json
import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import eigsh
import cell2 as C


def kron_list(ops):
    out = sps.identity(1, format='csr', dtype=complex)
    for o in ops:
        out = sps.kron(out, sps.csr_matrix(o), format='csr')
    return out


def run_operator(k, lam, par, with_left=True, with_right=True):
    """Factors: [B_left (4)] + [P (9)] * k + [B_right (4)]."""
    dims = [4] + [9] * k + [4]
    HL = C.h_pair('L', par['DN'], par['DC'], par['kap'], par['sg'], par['DCd'])     # (P, B)
    HR = C.h_pair('R', par['DN'], par['DC'], par['kap'], par['sg'], par['DCd'])     # (P, B)
    HRp = HR.reshape(9, 4, 9, 4).transpose(1, 0, 3, 2).reshape(36, 36)              # (B, P)
    V = np.einsum('acbd,ef,gh->aecgbfdh', (lam * C.port_V(par['Jp'], par['Jz'])).reshape(3, 3, 3, 3),
                  np.eye(3), np.eye(3)).reshape(81, 81)
    n = len(dims)
    H = sps.csr_matrix((int(np.prod(dims)),) * 2, dtype=complex)

    def two(op, i):          # op on factors (i, i+1)
        left = int(np.prod(dims[:i])); right = int(np.prod(dims[i + 2:]))
        return sps.kron(sps.kron(sps.identity(left), sps.csr_matrix(op)), sps.identity(right), format='csr')
    if with_left:
        H = H + two(HRp, 0)
    for i in range(1, k):
        H = H + two(V, i)
    if with_right:
        H = H + two(HL, k)
    return H


def lowest(H):
    if H.shape[0] <= 3000:
        return float(np.linalg.eigvalsh(H.toarray())[0])
    return float(eigsh(H, k=1, which='SA', tol=1e-10)[0][0])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--mu', type=float, default=1.0)
    ap.add_argument('--K', type=float, default=3.0)
    ap.add_argument('--lams', type=float, nargs='+', default=[0.0, 0.1, 0.2])
    ap.add_argument('--kmax', type=int, default=4)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    import effective2 as E
    par = dict(E.DEFAULT)
    vmin = float(np.linalg.eigvalsh(C.port_V(par['Jp'], par['Jz']))[0])
    print('frustrated shared port eps =', lowest(run_operator(1, 0.0, par)), ' min spec V =', vmin)
    res = {}
    for lam in a.lams:
        e = {k: lowest(run_operator(k, lam, par)) for k in range(1, a.kmax + 1)}
        f = lowest(run_operator(2, lam, par, with_right=False))
        u = {'0': -a.mu + a.K}
        u.update({str(k): (-a.mu + e[k]) / (k + 1) for k in e})
        u['tail'] = min((-a.mu + 2 * f + vmin * lam * (k - 3)) / (k + 1) for k in range(a.kmax + 1, 3000))
        m = {k: v - u['2'] for k, v in u.items() if k != '2'}
        print(f'lam={lam}: e_k={ {k: round(v, 5) for k, v in e.items()} }  f={f:.5f}')
        print('   u_k =', {k: round(v, 5) for k, v in u.items()}, ' min margin', round(min(m.values()), 5), 'at', min(m, key=m.get))
        res[str(lam)] = dict(e=e, f=f, u=u, margins=m)
    if a.out:
        json.dump(dict(args=vars(a), par={k: list(v) if isinstance(v, tuple) else v for k, v in par.items()}, vmin=vmin, res=res),
                  open(a.out, 'w'), indent=1)

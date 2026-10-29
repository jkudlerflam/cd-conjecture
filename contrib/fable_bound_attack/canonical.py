#!/usr/bin/env python3
"""Fixed-density Casimir coefficients for chains with conserved particle numbers, and the
bond scale M of the bond that makes that density the ground state.

For a family whose bond at zero chemical potentials is h_0, with filling charges N_a, we
compute E(T, L) canonically in the sector T = (n_a L) for commensurate L, fit A from those
energies (the chemical-potential term is then a constant and drops out of A), estimate
mu_a = [E(T + e_a) - E(T - e_a)]/2 at the largest L (so that h = h_0 - sum_a mu_a n_a has T
as its ground sector up to 1/L corrections), and report M(h).  Spin-like charge columns are
restricted to the two smallest |q|.

    python3 canonical.py --family tj --dens 0.6667 --Ls 6 9 12 --out data/tj_canon.jsonl
    python3 canonical.py --family gperm --par "par=[0,0,0,1]" sign=1 --dens 0.25 0.25 0.5 --Ls 4 8 12 ...
"""
import argparse, json, time
import numpy as np
import families as F
from symring import symring
from core import hnorm, fit_cas, LOG2
from scan import parse_par

# family -> (constructor at zero mu, filling columns, symmetric columns, number operators per filling column)
def setup(family, par):
    if family == 'tj':
        h, d, ch = F.tj(0.0, sign=par.get('sign', -1))
        nops = [np.diag([0.0, 1, 1])]
        return h, d, ch, [0], [1], nops, [0, 1, 1]
    if family == 'osp12':
        h, d, ch = F.osp12(par['theta'], 0.0)
        # charges: only the fermion S_z-like column; add a filling column (fermion number)
        ch = np.array([[0, 0], [1, 1], [1, -1]])
        return h, d, ch, [0], [1], [np.diag([0.0, 1, 1])], [0, 1, 1]
    if family == 'hubbard':
        h, d, ch = F.hubbard(par['U'], 0.0, par.get('t', 1.0))
        # columns (N_up, N_dn): use N = N_up + N_dn as filling, S_z = N_up - N_dn as symmetric
        ch = np.array([[0, 0], [1, 1], [1, -1], [2, 0]])
        n = np.diag([0.0, 1, 1, 2])
        return h, d, ch, [0], [1], [n], [0, 1, 1, 0]
    if family == 'gperm':
        parity = par['par']; sign = par.get('sign', 1)
        d = len(parity)
        h, d, ch = F.graded_perm(parity, [0.0] * d, sign)
        nops = []
        for a in range(1, d):
            n = np.zeros((d, d)); n[a, a] = 1; nops.append(n)
        return h, d, ch, list(range(d - 1)), [], nops, list(parity)
    raise ValueError(family)


def twisted(h, d, parity, L, phi):
    """(1 (x) R) h (1 (x) R^dag), R = exp(i phi n_F / L): a uniform gauge twist equivalent to a
    total flux phi through the ring for the fermionic species (parity 1)."""
    R = np.diag(np.exp(1j * phi * np.asarray(parity, float) / L))
    U = np.kron(np.eye(d), R)
    return U @ h @ U.conj().T


def best_twist(h0, d, parity, L, ch, targets, nlev, tol, entropy):
    best = None
    for phi in (0.0, np.pi):
        r = symring(twisted(h0, d, parity, L, phi), d, L, ch, targets=targets, nlev=nlev, tol=tol, entropy=entropy)
        r['phi'] = phi
        if best is None or r['E0'] < best['E0'] - 1e-12:
            best = r
    return best


def run(family, par, dens, Ls, nlev=2, tol=1e-10, entropy=True):
    h0, d, ch, fillcols, symcols, nops, parity = setup(family, par)
    ch = np.asarray(ch)
    t0 = time.time()
    rows = []
    for L in Ls:
        T = np.zeros(ch.shape[1], dtype=int)
        for c, n in zip(fillcols, dens):
            T[c] = int(round(n * L))
        targets = [tuple(T)]
        for c in symcols:
            for q in (1, -1):
                T2 = T.copy(); T2[c] = q; targets.append(tuple(T2))
        r = best_twist(h0, d, parity, L, ch, targets, nlev, tol, entropy)
        rows.append(r)
    E0 = [r['E0'] for r in rows]
    # chemical potentials at the largest L
    L = Ls[-1]
    T = np.array(rows[-1]['target'])
    mus = []
    for c in fillcols:
        Ep, Em = [], []
        for s, store in ((1, Ep), (-1, Em)):
            T2 = T.copy(); T2[c] += s
            targets = [tuple(T2)]
            for cc in symcols:
                for q in (1, -1, 0):
                    T3 = T2.copy(); T3[cc] = q; targets.append(tuple(T3))
            # parity of S_z must match particle number parity; include 0 and +-1
            rr = best_twist(h0, d, parity, L, ch, targets, 1, tol, False)
            store.append(rr['E0'])
        mus.append((Ep[0] - Em[0]) / 2)
    h = h0.copy()
    for mu, n in zip(mus, nops):
        h = h - mu * F.onsite(n, d)
    M = float(hnorm(h, d)); M0 = float(hnorm(h0, d))
    A, e0, b = fit_cas(Ls, E0)
    A2 = fit_cas(Ls[-2:], E0[-2:], quartic=False)[0]
    res = {'family': family, 'par': par, 'dens': list(dens), 'd': d, 'Ls': list(Ls), 'E0': E0,
           'k0': [r['k0'] for r in rows], 'phi': [r['phi'] for r in rows], 'target': [list(r['target']) for r in rows],
           'gapL': [r['gap'] * r['L'] for r in rows], 'v_est': [r['v_est'] for r in rows],
           'mu': mus, 'M': M, 'M_mu0': M0, 'A_fit': A, 'e0_fit': e0, 'A_last2': A2,
           'ratio': A / (np.pi * M * np.log(d) / LOG2), 'ratio_last2': A2 / (np.pi * M * np.log(d) / LOG2),
           'ratio_mu0': A / (np.pi * M0 * np.log(d) / LOG2), 'time': time.time() - t0}
    if len(Ls) >= 3:
        res['A_last3'] = fit_cas(Ls[-3:], E0[-3:])[0]
    if entropy and len(Ls) >= 2:
        S = [r['S_half'] for r in rows]
        res['S_half'] = S
        res['c_half'] = 3 * (S[-1] - S[-2]) / np.log(Ls[-1] / Ls[-2])
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--family', required=True)
    ap.add_argument('--par', nargs='*', default=[])
    ap.add_argument('--dens', nargs='+', type=float, required=True)
    ap.add_argument('--Ls', nargs='+', type=int, required=True)
    ap.add_argument('--nlev', type=int, default=2)
    ap.add_argument('--no-entropy', action='store_true')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    res = run(a.family, parse_par(a.par), a.dens, a.Ls, nlev=a.nlev, entropy=not a.no_entropy)
    with open(a.out, 'a') as f:
        f.write(json.dumps(res) + '\n')
    print(json.dumps({k: res[k] for k in ('family', 'par', 'dens', 'mu', 'M', 'A_fit', 'A_last2', 'ratio', 'ratio_last2', 'ratio_mu0', 'c_half', 'gapL', 'time') if k in res}))

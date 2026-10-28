#!/usr/bin/env python3
"""First-order spin-charge-separated spectrum of a vacancy ring, and comparison with exact ED (sqed.py).

Prediction (exact first-order degenerate perturbation theory in J about the free charge states):
  E = t E_chi(theta_r) + J alpha_chi eps_{r,i} - mu_v M,
where eps_{r,i} are the frozen-ring levels of the word on M squeezed sites with R-eigenvalue
exp(i theta_r), theta_r = 2 pi r / M (sqed.Ring with L = M), chi runs over free hard-core-boson states
with flux theta_r (fermion momenta k_a = (2 pi n_a + theta_r + pi (M - 1)) / L), alpha_chi =
(L/M)(rho^2 - |G|^2), G = (1/L) sum_a exp(i k_a), and the physical momentum is sum_a k_a.

usage: separated.py --L 10 --word BPPBPP --Q2 0 --J 0.25 --mu_v 0.6 [--nint 6 --nlev 6] --out f.json
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sqed  # noqa: E402


def charge_states(L, M, theta, emax, t=1.0):
    """Free hard-core-boson states (as fermions) with flux theta, energy within emax of the lowest."""
    base = []
    phi = theta + np.pi * (M - 1)
    ks = (2 * np.pi * np.arange(L) + phi) / L
    eps = -2 * t * np.cos(ks)
    order = np.argsort(eps)
    out = []
    # all M-subsets of the lowest min(L, M + 6) orbitals
    pool = order[:min(L, M + 6)]
    e0 = np.sort(eps)[:M].sum()
    for sub in itertools.combinations(pool, M):
        sub = np.array(sub)
        E = eps[sub].sum()
        if E > e0 + emax:
            continue
        kk = ks[sub]
        P = kk.sum()
        G = np.exp(1j * kk).sum() / L
        rho = M / L
        al = (L / M) * (rho ** 2 - abs(G) ** 2)
        mphys = int(round(P * L / (2 * np.pi))) % L
        out.append((float(E), min(mphys, L - mphys), float(al), float(E - e0)))
    return out


def predicted(L, word, Q2, t, J, mu_v, lam, nint=6, emax=1.0):
    M = len(word)
    froz = sqed.Ring(M, word, Q2, 0.0, 1.0, 0.0, lam)
    levels = {}
    for r in range(M // 2 + 1):
        pair = not (r == 0 or 2 * r == M)
        w = np.sort(froz.lowest(r, nev=(2 * nint if pair else nint), tol=1e-12))
        eps = w[::2] if pair else w                  # +-theta pairs appear twice in the cos sector
        thetas = [2 * np.pi * r / M] + ([-2 * np.pi * r / M] if pair else [])
        for th in thetas:
            for (E, mph, al, ech) in charge_states(L, M, th, emax, t):
                for i, e in enumerate(eps):
                    levels.setdefault(mph, []).append((E + J * al * e - mu_v * M, r, i, ech))
    return {m: sorted(v) for m, v in levels.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--word', required=True)
    ap.add_argument('--Q2', type=int, default=0)
    ap.add_argument('--t', type=float, default=1.0)
    ap.add_argument('--J', type=float, nargs='+', default=[0.25])
    ap.add_argument('--mu_v', type=float, default=0.6)
    ap.add_argument('--lam', type=float, default=0.15)
    ap.add_argument('--nlev', type=int, default=6)
    ap.add_argument('--out', default=None)
    ap.add_argument('--emax', type=float, default=1.0)
    a = ap.parse_args()
    L, M = a.L, len(a.word)
    res = dict(args=vars(a), runs=[])
    for J in a.J:
        pred = predicted(L, a.word, a.Q2, a.t, J, a.mu_v, a.lam, nint=a.nlev, emax=a.emax)
        ring = sqed.Ring(L, a.word, a.Q2, a.t, J, a.mu_v, a.lam)
        run = dict(J=J, sectors={})
        for m in range(L // 2 + 1):
            ed = ring.lowest(m, nev=a.nlev, tol=1e-12)
            lab = pred.get(m, [])[:len(ed)]
            pr = np.array([x[0] for x in lab] + [np.nan] * (len(ed) - len(lab)))
            run['sectors'][str(m)] = dict(ed=[float(x) for x in ed], pred=[float(x) for x in pr],
                                          labels=[dict(r=x[1], i=x[2], charge_exc=x[3]) for x in lab])
            print(f'J={J} L={L} {a.word} 2Q={a.Q2} m={m} ED={np.round(ed, 6).tolist()} '
                  f'pred={np.round(pr, 6).tolist()}', flush=True)
        res['runs'].append(run)
        if a.out:
            json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Tables for contrib/torus_routeB from the JSON outputs of run_vumps.py and run_obc.py.

  python3 analyze.py fes data/fes_*.json      finite-entanglement scaling
  python3 analyze.py obc data/obc_*.json      open-chain entropy fits
  python3 analyze.py exact                    exact free-fermion open-chain fits (calibration)

FES: c_FES = 6 dS/d ln xi from a linear fit over the last three chi (S averaged over the bonds
of the unit cell, xi in sites), and c_kappa from the exponent of xi ~ chi^kappa over the same chi.  OBC: least squares fit of
  S(l) = (c/6) ln[(2L/pi) sin(pi l/L)] + A + sum_a [B_a cos(2 k_a l) + C_a sin(2 k_a l)] / chord(l)
over 3 <= l <= L-3, with the free Fermi momenta k_a (one frequency for equal momenta).
Written by Claude (Opus 5.5) for contrib/torus_routeB, 2026-10-27.
"""
import sys, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def c_of_kappa(kap):
    """Invert kappa = 6 / (c (sqrt(12/c) + 1)) (Pollmann, Mukerjee, Turner, Moore 2009)."""
    from scipy.optimize import brentq
    f = lambda c: 6 / (c * (np.sqrt(12 / c) + 1)) - kap
    try:
        return brentq(f, 0.05, 60)
    except ValueError:
        return np.nan


def fes(files):
    print('c_FES = 6 dS/dln xi and c_kappa from xi ~ chi^kappa, both over the last three chi')
    print(f"{'file':40s} {'chi':>5s} {'S':>8s} {'xi':>8s} {'c_FES':>7s} {'c_kap':>7s} {'c_pair':>8s}  densities")
    for fn in files:
        d = json.load(open(fn))
        R = d['results']
        if not R:
            continue
        S = np.array([np.mean(r['S']) for r in R])
        xi = np.array([r['xi_sites'] if 'xi_sites' in r else r['xi'] for r in R])
        chis = [r['chi'] for r in R]
        cs = [np.nan]
        for i in range(1, len(R)):
            cs.append(6 * (S[i] - S[i - 1]) / np.log(xi[i] / xi[i - 1]) if xi[i] != xi[i - 1] else np.nan)
        c3, ck = np.nan, np.nan
        if len(R) >= 3:
            p = np.polyfit(np.log(xi[-3:]), S[-3:], 1)
            c3 = 6 * p[0]
            kap = np.polyfit(np.log(np.array(chis[-3:], float)), np.log(xi[-3:]), 1)[0]
            ck = c_of_kappa(kap)
        for i, r in enumerate(R):
            n = np.round(np.mean(r['n'], 0), 3).tolist() if 'n' in r else ''
            last = i == len(R) - 1
            print(f"{os.path.basename(fn) if i == 0 else '':40s} {chis[i]:5d} {S[i]:8.4f} {xi[i]:8.2f} "
                  f"{(c3 if last else np.nan):7.2f} {(ck if last else np.nan):7.2f} {cs[i]:8.3f}  {n}")


def cc_fit(S, L, k, lmin=3):
    l = np.arange(1, L)
    sel = (l >= lmin) & (l <= L - lmin)
    chord = (2 * L / np.pi) * np.sin(np.pi * l / L)
    cols = [np.log(chord) / 6, np.ones_like(chord)]
    for kk in np.unique(np.round(k, 10)):
        cols += [np.cos(2 * kk * l) / chord, np.sin(2 * kk * l) / chord]
    A = np.array(cols).T[sel]
    coef, res, *_ = np.linalg.lstsq(A, np.asarray(S)[sel], rcond=None)
    return coef[0], float(np.sqrt(np.mean((A @ coef - np.asarray(S)[sel]) ** 2)))


def k_of(args):
    import tcases as M
    if args.get('mode') == 'model':
        return M.K_EQ if args['case'] in ('freeeq', 'iii') else M.K_DIST
    return None


def obc(files):
    print(f"{'file':40s} {'chi':>5s} {'trunc':>8s} {'S(L/2)':>8s} {'c_fit':>7s} {'rms':>7s}  densities")
    for fn in files:
        d = json.load(open(fn))
        a = d['args']
        L = a['L']
        k = k_of(a)
        for i, r in enumerate(d['recs']):
            S = r['S']
            if k is None:
                n = np.mean(r['n'], 0)
                k = np.pi * np.array(n)          # effective chains: momenta from densities
            c, rms = cc_fit(S, L, k)
            print(f"{os.path.basename(fn) if i == 0 else '':40s} {r['chi']:5d} {r['trunc']:8.1e} {S[L // 2 - 1]:8.4f} "
                  f"{c:7.3f} {rms:7.4f}  {np.round(np.mean(r['n'], 0), 3).tolist()}")


def exact():
    import tcases as M, fermions as F
    for name, k in (('distinct', M.K_DIST), ('equal', M.K_EQ)):
        t, mu = F.free_params(k)
        for L in (24, 32, 64, 128):
            S = np.zeros(L - 1)
            for a in range(len(k)):
                H = -t[a] * (np.eye(L, k=1) + np.eye(L, k=-1)) - mu[a] * np.eye(L)
                e, U = np.linalg.eigh(H)
                occ = U[:, e < 0]
                Cm = occ @ occ.conj().T
                for ell in range(1, L):
                    nu = np.linalg.eigvalsh(Cm[:ell, :ell]).clip(1e-14, 1 - 1e-14)
                    S[ell - 1] += -np.sum(nu * np.log(nu) + (1 - nu) * np.log(1 - nu))
            c, rms = cc_fit(S, L, k)
            print(f'exact free, {name} momenta, L = {L}: c_fit = {c:.3f} (rms {rms:.4f}), S(L/2) = {S[L // 2 - 1]:.4f}')


if __name__ == '__main__':
    if sys.argv[1] == 'fes':
        fes(sys.argv[2:])
    elif sys.argv[1] == 'obc':
        obc(sys.argv[2:])
    else:
        exact()

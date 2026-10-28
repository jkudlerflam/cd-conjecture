#!/usr/bin/env python3
"""Ground-state search over all competitive sectors of one ring, for the vacancy chain or the frozen chain.

Sectors are (M occupied, cyclic type word, 2Q, momentum m).  Candidate (M, word) pairs are screened with the
first-order energy E1 = t E_free(L, M) + J alpha_{L,M} E_run(word) - mu_v M, where E_free is the hard-core
boson ground energy, alpha_{L,M} the exact adjacency factor of the charge sea
(../counterexample_campaign/doping_check/README.md) and E_run = sum_B(-mu) + sum_runs e_k the run-energy
estimate of ../routeA_demonstrator/selector.py (a rigorous lower bound for the frozen word, tight up to the
O(lam^2) logical energy).  Words contain runs of k = 1..4 ports between consecutive B's (BB costs K = 2).
Every candidate within `window` of the best E1 (and at least the best `nmin`) is then diagonalized exactly
in Q2 in {M mod 2, M mod 2 + 2} and every momentum.  --frozen sets M = L, J = 1, no vacancies.
Output: all sector levels, the global ground state, its degeneracy (+-k pairs and +-Q counted), and gaps.
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'routeA_demonstrator'))
import sqed  # noqa: E402
import selector  # noqa: E402
import vacbond as VB  # noqa: E402


def necklaces(M, kmax=4):
    """Cyclic words of length M made of units B P^k, 1 <= k <= kmax, up to rotation."""
    out = set()

    def rec(prefix, left):
        if left == 0:
            out.add(sqed.canonical(prefix))
            return
        for k in range(1, kmax + 1):
            if k + 1 <= left:
                rec(prefix + 'B' + 'P' * k, left - k - 1)
    rec('', M)
    out.add('P' * M)            # no B at all: the port-only XXZ chain
    return sorted(out)


def alpha(L, M):
    if M >= L:
        return 1.0
    return M / L - np.sin(np.pi * M / L) ** 2 / (M * L * np.sin(np.pi / L) ** 2)


def run_energy(word, ek, mu, lam=0.15):
    """sum_B(-mu) + sum_runs e_k for a cyclic word with no BB; -lam M (a lower bound, every port-port
    bond has spectrum lam{-1,-1,-1,3}) for the word without B."""
    if 'B' not in word:
        return -lam * len(word)
    w = word
    i0 = w.index('B')
    w = w[i0:] + w[:i0]
    runs = [len(r) for r in w.split('B')[1:]]
    return -mu * w.count('B') + sum(ek[k] for k in runs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--t', type=float, default=1.0)
    ap.add_argument('--J', type=float, default=0.25)
    ap.add_argument('--mu_v', type=float, default=0.0)
    ap.add_argument('--lam', type=float, default=0.15)
    ap.add_argument('--frozen', action='store_true')
    ap.add_argument('--window', type=float, default=0.03)
    ap.add_argument('--nmin', type=int, default=3)
    ap.add_argument('--nev', type=int, default=4)
    ap.add_argument('--Mrange', type=int, nargs=2, default=None)
    ap.add_argument('--q0only', action='store_true', help='only the lowest |2Q| sector (M mod 2)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    L = a.L
    ek = selector.unit_energies(a.lam, VB.DEMO['D0'], VB.DEMO['D1'], VB.DEMO['kap'], VB.DEMO['Jp'], VB.DEMO['Jz'], kmax=4)
    mu = VB.DEMO['mu']
    cands = []
    Ms = [L] if a.frozen else range(a.Mrange[0] if a.Mrange else 2, (a.Mrange[1] if a.Mrange else L - 1) + 1)
    for M in Ms:
        if M < 2:
            continue
        ef = 0.0 if a.frozen else a.t * sqed.lowest_free(L, M, 0.0)
        J = 1.0 if a.frozen else a.J
        al = alpha(L, M)
        for w in necklaces(M):
            E1 = ef + J * al * run_energy(w, ek, mu, a.lam) - (0 if a.frozen else a.mu_v) * M
            cands.append(dict(M=M, word=w, E1=float(E1)))
    cands.sort(key=lambda c: c['E1'])
    e1min = cands[0]['E1']
    chosen = [c for i, c in enumerate(cands) if i < a.nmin or c['E1'] <= e1min + a.window]
    res = dict(args=vars(a), candidates=cands[:40], chosen=chosen, sectors=[])
    t0 = time.time()
    for c in chosen:
        M, w = c['M'], c['word']
        q0 = M % 2
        for Q2 in ((q0,) if a.q0only else (q0, q0 + 2)):
            if a.frozen:
                ring = sqed.Ring(L, w, Q2, 0.0, 1.0, 0.0, a.lam)
            else:
                ring = sqed.Ring(L, w, Q2, a.t, a.J, a.mu_v, a.lam)
            for m in range(L // 2 + 1):
                lev = ring.lowest(m, nev=a.nev, tol=1e-12)
                res['sectors'].append(dict(M=M, word=w, Q2=Q2, m=m, dim=[ring.D, ring.C], levels=[float(x) for x in lev]))
                print(f'L={L} M={M} {w} 2Q={Q2} m={m} dim={ring.D}x{ring.C} {np.round(lev, 9).tolist()} '
                      f't={time.time() - t0:.0f}s', flush=True)
            json.dump(res, open(a.out + '.tmp', 'w'), indent=1)
            os.replace(a.out + '.tmp', a.out)
    # ground state and degeneracy
    allv = []
    for s in res['sectors']:
        mult = 2 if s['Q2'] != 0 else 1
        for x in s['levels']:
            allv.append((x, mult, s['M'], s['word'], s['Q2'], s['m']))
    allv.sort()
    E0 = allv[0][0]
    tol = 1e-8 * max(1.0, abs(E0))
    deg = sum(v[1] for v in allv if v[0] < E0 + tol)
    gs = [dict(M=v[2], word=v[3], Q2=v[4], m=v[5]) for v in allv if v[0] < E0 + tol]
    higher = sorted(set(round(v[0] - E0, 12) for v in allv if v[0] >= E0 + tol))[:8]
    sat = [dict(M=s['M'], word=s['word'], Q2=s['Q2'], m=s['m']) for s in res['sectors']
           if max(s['levels']) < E0 + tol]
    res['summary'] = dict(E0=E0, degeneracy=deg, gs_sectors=gs, gaps=higher, saturated_sectors=sat,
                          time_s=time.time() - t0)
    print(json.dumps(res['summary']), flush=True)
    json.dump(res, open(a.out + '.tmp', 'w'), indent=1)
    os.replace(a.out + '.tmp', a.out)


if __name__ == '__main__':
    main()

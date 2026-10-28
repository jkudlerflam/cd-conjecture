#!/usr/bin/env python3
"""Ring estimators for a chain against its own free (Hartree-Fock) part, both on the same rings.
For every L: ground-state sector N, energy E0, ring-entanglement fit c_S(L) (S against
(1/3) ln[(L/pi) sin(pi l/L)]), and the lowest gaps in the ground sector.  The free part is solved
exactly (Gaussian) on the same ring with the same twist.  Casimir estimator: the difference
dE(L) = E0_chain(L) - E0_free(L) is fitted to A L + B / L, which removes the extensive part and the
shell (incommensurability) oscillations common to both, and c_Cas = n - 12 B / (pi sum_p |v_p|) with
the Hartree-Fock Fermi velocities.  usage: ring_fit.py RING.json SPEC6 [Lmin]
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28."""
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import models34 as Mo
import analyze34 as A
from gauss_calib import ring_free, fit_cS


def free_levels(T, M, L, N, const, k=6, m=8):
    th = 0.0 if (N - 1) % 2 == 0 else np.pi
    ks = (2 * np.pi * np.arange(L) + th) / L
    eps = np.sort(np.concatenate([np.linalg.eigvalsh(A.hk(T, M, q)) for q in ks]))
    occ = list(range(N)); emp = list(range(N, len(eps)))
    E0 = eps[:N].sum()
    top = occ[-m:]; bot = emp[:m]
    Es = [E0]
    for r in (1, 2, 3):
        for hs in itertools.combinations(top, r):
            for ps in itertools.combinations(bot, r):
                Es.append(E0 - eps[list(hs)].sum() + eps[list(ps)].sum())
    Es = np.sort(Es)[:k]
    return Es + L * const


def main():
    ring = json.load(open(sys.argv[1]))
    spec6 = sys.argv[2]
    Lmin = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    lev = ring['info']['level']
    n = ring['info']['n']
    if lev.startswith('p'):
        _, T, M, _ = Mo.proxy(spec6)
    else:
        _, _, T, M, _ = Mo.base(spec6)
    Kq = Mo.quad_bond(T, M, n)
    const = -float(np.trace(Kq).real) / len(Kq)
    vs = sum(abs(f['v']) for f in A.fermi(T, M))
    byL = {}
    for r in ring['results']:
        byL.setdefault(r['L'], []).append(r)
    rows = []
    for L in sorted(byL):
        rs = byL[L]
        g = min(rs, key=lambda r: r['levels'][0])
        Nf, Ef, Sf, EN = ring_free(T, M, L, const)
        fl = free_levels(T, M, L, g['N'], const)
        gaps_c = np.array(g['levels'][1:4]) - g['levels'][0]
        gaps_f = fl[1:4] - fl[0]
        rows.append(dict(L=L, N=g['N'], N_free=Nf, E0=g['levels'][0], E0_free=Ef, E0_free_sameN=EN[g['N']],
                         cS=fit_cS(g['S'], L), cS_free=fit_cS(Sf, L) if Nf == g['N'] else None,
                         gaps=gaps_c.tolist(), gaps_free_sameN=gaps_f.tolist()))
    Ls = np.array([r['L'] for r in rows if r['L'] >= Lmin], float)
    dE = np.array([r['E0'] - r['E0_free_sameN'] for r in rows if r['L'] >= Lmin])
    Am = np.vstack([Ls, 1 / Ls]).T
    (a, b), *_ = np.linalg.lstsq(Am, dE, rcond=None)
    cCas = n - 12 * b / (np.pi * vs)
    out = dict(spec=ring['info']['spec'], rows=rows, v_sum=vs, fit_A=float(a), fit_B=float(b), c_Cas=float(cCas))
    for r in rows:
        print('L=%d N=%d (free %d) cS=%.3f cS_free=%s  gaps %s  free %s' % (
            r['L'], r['N'], r['N_free'], r['cS'], None if r['cS_free'] is None else round(r['cS_free'], 3),
            np.round(r['gaps'], 4).tolist(), np.round(r['gaps_free_sameN'], 4).tolist()))
    print('Casimir difference fit over L >= %d: A=%.4e B=%.4e  c_Cas = %.3f' % (Lmin, a, b, cCas))
    json.dump(out, open(sys.argv[1].replace('.json', '_fit.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()

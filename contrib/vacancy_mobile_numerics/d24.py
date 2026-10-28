#!/usr/bin/env python3
"""The 23-letter route-A bond at n = 1 plus one vacancy (d = 24): small-ring ground states.

Bond: ../explicit_candidate/bond23.py (literal 529 x 529, complex), lifted exactly as in
../counterexample_campaign/doping/bond24.py (vacancy hop -t_v, J h on occupied pairs, -mu_v per occupied
site).  Conserved and used: occupied number, the cyclic word of site types D (data letters) / M (marker
letters), the charge Q, the arrow-legality of every squeezed bond (moves preserve outer arrows, so the
legal sector is an exact symmetry sector), and momentum (complex projector; the bond is complex).

usage: d24.py check                     literal 576 x 576 lift vs squeezed ED on small rings
       d24.py frozen --L 6              frozen d = 23 ring (no vacancies)
       d24.py vac --L 7 --J 0.25 --mu_v 0.0
       options: --point fast|working  (coefficient point at n = 1, see POINTS)
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, time, argparse, itertools, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'explicit_candidate'))
import sqed  # noqa: E402
import bond23 as B  # noqa: E402

# Two coefficient points at n = 1 (cell length b = 2, mu_sel = 2t - 1.9 = 0.1, the selector midpoint).
#  working: the checked working point of ../explicit_candidate/README.md Section 1.3 / 2.6 (D_x = 4, g = 0.05),
#           whose logical hopping is |T| ~ 4e-7 lambda^2: the payload is frozen at any resolvable scale.
#  fast:    same terms, coefficients chosen so the logical hopping is ~1e-3 lambda^2 (D_x = D_M = 0.5,
#           g = 0.3, Delta = 0.4 >= g); the code is still the cell ground space (gap 0.2).
POINTS = dict(
    working=dict(gp=0.05, gm=0.05, Dx=4.0, DM=4.0, Delta=0.0, td=1e-3, musel=0.1, K=10.0,
                 phx={'0': 0.4, '+': 0.0, '-': 0.0}, phd={'0': 0.9, '+': 0.0, '-': 0.0}),
    fast=dict(gp=0.3, gm=0.3, Dx=0.5, DM=0.5, Delta=0.4, td=1e-3, musel=0.1, K=10.0,
              phx={'0': np.pi / 2, '+': 0.0, '-': 0.0}, phd={'0': 0.9, '+': 0.0, '-': 0.0}))
LAM = 0.2


def params(point, lam=LAM, eps=None):
    P = B.params(lam=lam, **POINTS[point])
    if eps is None:
        T, _, _ = B.onebody_T(1, P)
        eps = 0.5 * abs(T[0, 0]) * lam ** 2         # data field = lambda^2 * (1/2) lambda_min(R_s), README rule
    P['eps'] = eps
    return P


def alphabet(P):
    h = B.bond(P)
    D = B.D
    marker = np.array([l[0] in ('Me', 'M') for l in B.LET])
    letters = {'D': np.where(~marker)[0], 'M': np.where(marker)[0]}
    tp = marker.astype(int)
    T2 = np.add.outer(2 * tp, tp).reshape(-1)
    assert np.abs(h * (T2[:, None] != T2[None, :])).max() == 0, 'marker type not conserved'
    hb = {}
    for a in 'DM':
        for b in 'DM':
            idx = np.array([x * D + y for x in letters[a] for y in letters[b]])
            hb[a + b] = h[np.ix_(idx, idx)]
    q = np.array([B.charge(l) for l in B.LET], int)
    allowed = np.array([[B.ARROW[l[0]][1] == B.ARROW[r[0]][0] for r in B.LET] for l in B.LET])
    return dict(letters=letters, q2=q, hb=hb, allowed=allowed), h


def ring(L, word, Q, t, J, mu_v, P, alph=None):
    if alph is None:
        alph, _ = alphabet(P)
    return sqed.Ring(L, word, Q, t, J, mu_v, P['lam'], alph=alph)


def spectrum(r, nev=3):
    out = {}
    for m in range(r.L):
        out[m] = [float(x) for x in r.lowest(m, nev=nev, tol=1e-12)]
    return out


def check(P):
    """Literal lift on the full physical basis (legal squeezed bonds only) vs squeezed ED."""
    spec = importlib.util.spec_from_file_location(
        'bond24', os.path.join(HERE, '..', 'counterexample_campaign', 'doping', 'bond24.py'))
    b24 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b24)
    alph, h = alphabet(P)
    out = []
    for (L, word, Q, t, J, mu_v) in [(5, 'DMDM', 1, 1.0, 0.4, 0.1), (6, 'DMDM', 2, 1.0, 0.3, 0.0),
                                     (5, 'DDM', 1, 1.0, 0.5, 0.2)]:
        h24 = b24.vacancy_lift(h, hopping=t, internal_scale=J, chemical_potential=mu_v)
        Mo = len(word)
        rots = sorted(set(word[i:] + word[:i] for i in range(Mo)))
        states = []
        for x in itertools.combinations(range(L), Mo):
            for w in rots:
                for lab in itertools.product(*[alph['letters'][ty] for ty in w]):
                    if sum(alph['q2'][a] for a in lab) != Q:
                        continue
                    if not all(alph['allowed'][lab[j], lab[(j + 1) % Mo]] for j in range(Mo)):
                        continue
                    s = [0] * L
                    for j, xi in enumerate(x):
                        s[xi] = 1 + lab[j]
                    states.append(tuple(s))
        idx = {s: i for i, s in enumerate(states)}
        n = len(states)
        H = np.zeros((n, n), complex)
        d = 24
        for i, s in enumerate(states):
            for b in range(L):
                b1 = (b + 1) % L
                k = s[b] * d + s[b1]
                for o in np.nonzero(h24[:, k])[0]:
                    ns = list(s); ns[b], ns[b1] = divmod(o, d)
                    H[idx[tuple(ns)], i] += h24[o, k]        # KeyError = leak out of the sector
        wb = np.linalg.eigvalsh(H)[:8]
        r = sqed.Ring(L, word, Q, t, J, mu_v, P['lam'], alph=alph)
        ws = np.sort(np.concatenate([r.lowest(m, nev=8, tol=1e-12) for m in range(L)]))[:8]
        out.append(dict(L=L, word=word, Q=Q, dim=n, dim_sqed=r.N, herm=float(np.abs(H - H.conj().T).max()),
                        max_diff=float(np.abs(wb - ws).max())))
        print(out[-1], flush=True)
    return out


def words(M, kmax=3):
    """Cyclic D/M words made of cells D^k M, 1 <= k <= kmax."""
    out = set()

    def rec(prefix, left):
        if left == 0:
            out.add(sqed.canonical(prefix)); return
        for k in range(1, kmax + 1):
            if k + 1 <= left:
                rec(prefix + 'D' * k + 'M', left - k - 1)
    rec('', M)
    return sorted(out)


def frozen_energy(word, P, alph, Qs):
    best = np.inf
    for Q in Qs:
        r = sqed.Ring(len(word), word, Q, 0.0, 1.0, 0.0, P['lam'], alph=alph)
        if r.N == 0:
            continue
        for m in range(len(word)):
            best = min(best, r.lowest(m, nev=1, tol=1e-12)[0])
    return best


def scan(L, P, frozen, t, J, mu_v, nmin=3, window=0.02, nev=4):
    alph, _ = alphabet(P)
    Ms = [L] if frozen else range(2, L)
    cands = []
    eint = {}
    for M in Ms:
        for w in words(M):
            if w not in eint:
                ncell = w.count('M')
                eint[w] = frozen_energy(w, P, alph, range(0, ncell + 2))
            if frozen:
                E1 = eint[w]
            else:
                al = M / L - np.sin(np.pi * M / L) ** 2 / (M * L * np.sin(np.pi / L) ** 2)
                E1 = t * sqed.lowest_free(L, M, 0.0) + J * al * eint[w] - mu_v * M
            cands.append(dict(M=M, word=w, E1=float(E1)))
    cands.sort(key=lambda c: c['E1'])
    chosen = [c for i, c in enumerate(cands) if i < nmin or c['E1'] <= cands[0]['E1'] + window]
    sectors = []
    for c in chosen:
        ncell = c['word'].count('M')
        for Q in range(0, ncell + 2):
            r = sqed.Ring(L, c['word'], Q, 0.0 if frozen else t, 1.0 if frozen else J, 0.0 if frozen else mu_v,
                          P['lam'], alph=alph)
            if r.N == 0:
                continue
            for m in range(L):
                lev = r.lowest(m, nev=nev, tol=1e-12)
                sectors.append(dict(M=c['M'], word=c['word'], Q=Q, m=m, dim=[r.D, r.C], levels=[float(x) for x in lev]))
            print(f"L={L} {c['word']} Q={Q} done, lowest {min(s['levels'][0] for s in sectors if s['word'] == c['word'] and s['Q'] == Q):.10f}",
                  flush=True)
    allv = sorted((x, s['M'], s['word'], s['Q'], s['m']) for s in sectors for x in s['levels'])
    E0 = allv[0][0]
    tol = 1e-9 * max(1, abs(E0))
    gs = [dict(M=v[1], word=v[2], Q=v[3], m=v[4]) for v in allv if v[0] < E0 + tol]
    gaps = sorted(set(round(v[0] - E0, 12) for v in allv if v[0] >= E0 + tol))[:8]
    return dict(candidates=cands[:20], chosen=chosen, sectors=sectors,
                summary=dict(E0=E0, degeneracy=len(gs), gs=gs, gaps=gaps))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['check', 'frozen', 'vac', 'info'])
    ap.add_argument('--L', type=int, default=6)
    ap.add_argument('--point', default='fast')
    ap.add_argument('--t', type=float, default=1.0)
    ap.add_argument('--J', type=float, default=0.25)
    ap.add_argument('--mu_v', type=float, default=0.0)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    P = params(a.point)
    t0 = time.time()
    if a.mode == 'info':
        T, E0, res = B.onebody_T(1, P)
        c = B.cell_from_bond(1, P)
        w = np.linalg.eigvalsh(c['H'])
        res = dict(point=a.point, T=[float(T[0, 0].real), float(T[0, 0].imag)], eps=P['eps'], lam=P['lam'],
                   cell_levels=[float(x) for x in w[:6]], checks=B.bond_checks(P))
    elif a.mode == 'check':
        res = check(P)
    else:
        res = scan(a.L, P, a.mode == 'frozen', a.t, a.J, a.mu_v)
        res['args'] = vars(a)
    res = dict(result=res, time_s=time.time() - t0, point=a.point,
               params={k: (v if not isinstance(v, dict) else dict(v)) for k, v in P.items()})
    s = json.dumps(res, indent=1, default=float)
    if a.out:
        open(a.out, 'w').write(s)
    print(json.dumps(res['result'].get('summary', res['result']) if isinstance(res['result'], dict) else res['result'],
                     default=float)[:2000])


if __name__ == '__main__':
    main()

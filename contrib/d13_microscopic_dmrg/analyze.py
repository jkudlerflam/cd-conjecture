#!/usr/bin/env python3
"""Calibrated central charges from data/obc/*.json.  Claude (Opus 5.5) for Jonah Kudler-Flam, 2026-10-29.

Estimator A (difference, as in ../vacancy_n2_dmrg): c = m + slope of S_model(l) - S_ref(l), where S_ref
is the exact profile of m independent open free-fermion chains with the sector's particle numbers.
  mobile member, n flavours: m = n + 1 chains of length L with (N, Nlog) particles [vacancy fluid + payload]
  frozen control (N = L): sites are cells, S sampled at cell boundaries l = b j, reference m = n chains
  of Mk sites with Nlog particles, chord variable in cells.
Estimator B (ratio): m c_fit(model) / c_fit(ref) with the same fit function and window.
Estimator C (half-chain pairs, mobile only): 6 [S(L2/2) - S(L1/2)] / ln[(L2+1)/(L1+1)], times m / (same for ref).
Writes data/summary.json.
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from fitlib import profile_fit, chord


def free_profile(L, N):
    """Exact S(l), l = 0..L, of N free fermions on an open chain of L sites (ground state)."""
    S = np.zeros(L + 1)
    if N == 0 or N == L:
        return S
    j = np.arange(1, L + 1)
    k = np.arange(1, N + 1)
    phi = np.sqrt(2 / (L + 1)) * np.sin(np.pi * np.outer(j, k) / (L + 1))
    C = phi @ phi.T
    for l in range(1, L):
        x = np.linalg.eigvalsh(C[:l, :l]); x = x[(x > 1e-13) & (x < 1 - 1e-13)]
        S[l] = -np.sum(x * np.log(x) + (1 - x) * np.log(1 - x))
    return S


def fit_diff(Sm, Sr, L, qs, m, lmin):
    l = np.arange(lmin, L - lmin + 1)
    cd, rms = profile_fit(np.asarray(Sm) - np.asarray(Sr), L, qs, lmin=lmin)
    return m + cd


def analyse(d):
    a = d['args']; L = a['L']; N = a['N']; Nl = a['Nlog']; p = a['p']; b = p + 1; Mk = d['Mk']
    n = p - 1
    rec = dict(name=d['name'], L=L, N=N, Nlog=Nl, p=p, Mk=Mk, rows=[])
    for r in d['results']:
        S = np.array(r['S'])
        dens = np.array(r['dens'])          # per site <Q>, <occ>, <marker>
        row = dict(chi=r['chi'], maxchi=r['maxchi'], trunc=r['trunc'], E=r['E'], S_mid=float(S[L // 2]),
                   vac_density=float(1 - dens[:, 1].mean()), marker_per_occ=float(dens[:, 2].sum() / max(dens[:, 1].sum(), 1e-12)))
        if N == L:     # frozen control: cell boundaries
            Sc = S[b * np.arange(Mk + 1)]
            ref = free_profile(Mk, Nl) * n if n == 1 else None
            if n == 1:
                qs = [2 * np.pi * Nl / (Mk + 1)]
                cA = []; cB = []
                for lmin in (1, 2):
                    if Mk - 2 * lmin < 3: continue
                    cd, _ = profile_fit(Sc - ref, Mk, qs, lmin=lmin); cA.append(1 + cd)
                    cm, _ = profile_fit(Sc, Mk, qs, lmin=lmin); cr, _ = profile_fit(ref, Mk, qs, lmin=lmin)
                    cB.append(cm / cr)
                row.update(cA=cA, cB=cB)
            else:
                cm, _ = profile_fit(Sc, Mk, [2 * np.pi * Nl / (Mk + 1)], lmin=1)
                row.update(c_raw=cm)
        else:
            m = n + 1
            ref = free_profile(L, N) + n * free_profile(L, Nl)  # payload modes share the logical count (n = 1)
            qs = [2 * np.pi * N / (L + 1), 2 * np.pi * Nl / (L + 1), 2 * np.pi * Mk / (L + 1)]
            cA = []; cB = []
            for lmin in (4, 6, 8):
                if L - 2 * lmin < 10: continue
                cd, _ = profile_fit(S - ref, L, qs, lmin=lmin); cA.append(m + cd)
                cm, _ = profile_fit(S, L, qs, lmin=lmin); cr, _ = profile_fit(ref, L, qs, lmin=lmin)
                cB.append(m * cm / cr)
            row.update(cA=cA, cB=cB, S_ref_mid=float(ref[L // 2]))
        rec['rows'].append(row)
    return rec


def main():
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, 'data', 'obc', '*.json'))):
        d = json.load(open(f))
        if not d['results']:
            continue
        r = analyse(d); out.append(r)
        print(r['name'], 'L', r['L'], 'N', r['N'], 'Nlog', r['Nlog'], 'Mk', r['Mk'])
        for w in r['rows']:
            extra = ' '.join(f'{k}={np.round(v, 3).tolist()}' for k, v in w.items() if k in ('cA', 'cB', 'c_raw'))
            print(f"  chi={w['chi']} maxchi={w['maxchi']} trunc={w['trunc']:.1e} E={w['E']:.10f} Smid={w['S_mid']:.4f} "
                  f"vac={w['vac_density']:.4f} mk/occ={w['marker_per_occ']:.4f} {extra}")
    # estimator C: half-chain pairs of mobile n = 1 runs at the same J
    mob = {r['L']: r for r in out if r['p'] == 2 and r['N'] < r['L'] and r['name'].startswith('n1_')}
    Ls = sorted(mob)
    pairs = []
    for L1, L2 in zip(Ls, Ls[1:]):
        a1, a2 = mob[L1]['rows'][-1], mob[L2]['rows'][-1]
        r1 = free_profile(L1, mob[L1]['N'])[L1 // 2] + free_profile(L1, mob[L1]['Nlog'])[L1 // 2]
        r2 = free_profile(L2, mob[L2]['N'])[L2 // 2] + free_profile(L2, mob[L2]['Nlog'])[L2 // 2]
        den = np.log((L2 + 1) / (L1 + 1))
        cm = 6 * (a2['S_mid'] - a1['S_mid']) / den; cr = 6 * (r2 - r1) / den
        pairs.append(dict(L1=L1, L2=L2, chi1=a1['chi'], chi2=a2['chi'], c_raw=cm, c_ref=cr, c_cal=2 * cm / cr))
        print('pair', pairs[-1])
    json.dump(dict(runs=out, pairs=pairs), open(os.path.join(HERE, 'data', 'summary.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()

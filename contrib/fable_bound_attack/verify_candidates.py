#!/usr/bin/env python3
"""Recompute the best optimizer candidates on larger rings with criticality diagnostics.

    python3 verify_candidates.py data/rescore_top.json --d 2 --n 3 --Ls 8 10 12 14 16 18 --out data/verify_d2.json

For each candidate bond: E_0(L), the two lowest gaps times L, the half-ring entropy and c_half,
the full Casimir fit, the two-point estimates and the drift score, M, and the spectrum of the
bond.  The ground-state momentum is read off by applying the translation operator."""
import argparse, json
import numpy as np
from core import Ring, hnorm, fit_cas, half_entropy, LOG2, spread


def translate(psi, d, L):
    return np.moveaxis(psi.reshape([d] * L), 0, -1).reshape(-1)


def analyse(h, d, Ls, nlev=3):
    out = {'Ls': list(Ls), 'E0': [], 'gapL1': [], 'gapL2': [], 'S_half': [], 'k0': []}
    for L in Ls:
        r = Ring(h, d, L)
        w, v = r.lowest(k=nlev, tol=1e-9, return_vec=True)
        out['E0'].append(float(w[0]))
        out['gapL1'].append(float((w[1] - w[0]) * L))
        out['gapL2'].append(float((w[2] - w[0]) * L) if len(w) > 2 else None)
        out['S_half'].append(half_entropy(v, d, L))
        tv = translate(v, d, L)
        ph = np.vdot(v, tv)
        out['k0'].append(float(np.angle(ph) * L / (2 * np.pi)) if abs(ph) > 0.5 else None)
    M = float(hnorm(h, d))
    den = np.pi * M * np.log(d) / LOG2
    A = fit_cas(Ls, out['E0'])[0]
    two = [fit_cas(Ls[i:i + 2], out['E0'][i:i + 2], quartic=False)[0] / den for i in range(len(Ls) - 1)]
    out.update({'M': M, 'A_fit': float(A), 'ratio_fit': float(A / den), 'two_point': [float(t) for t in two],
                'score': float(min(two) - 0.5 * (max(two) - min(two))),
                'c_half_pairs': [float(3 * (out['S_half'][i + 1] - out['S_half'][i]) / np.log(Ls[i + 1] / Ls[i])) for i in range(len(Ls) - 1)],
                'bond_spectrum': [float(x) for x in np.linalg.eigvalsh(h)], 'half_spread': spread(h) / 2})
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('top')
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--n', type=int, default=3)
    ap.add_argument('--Ls', nargs='+', type=int, required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    top = json.load(open(a.top))[str(a.d)]
    # distinct candidates (by rounded h or x)
    seen, cands = set(), []
    for r in top:
        if r.get('h_re') is None:
            continue
        key = tuple(np.round(np.array(r['h_re']).ravel(), 3))
        if key in seen:
            continue
        seen.add(key); cands.append(r)
        if len(cands) >= a.n:
            break
    res = []
    for r in cands:
        h = np.array(r['h_re']) + 1j * np.array(r['h_im'] if r.get('h_im') else np.zeros_like(r['h_re']))
        if np.allclose(h.imag, 0):
            h = h.real
        an = analyse(h, a.d, a.Ls)
        an['source'] = {k: r[k] for k in ('file', 'start', 'init', 'score_new')}
        res.append(an)
        print(json.dumps({k: an[k] for k in ('source', 'M', 'ratio_fit', 'two_point', 'score', 'gapL1', 'c_half_pairs', 'k0')}))
    json.dump(res, open(a.out, 'w'), indent=1)

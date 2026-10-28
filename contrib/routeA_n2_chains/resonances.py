#!/usr/bin/env python3
"""Fermi points, velocities and resonance margins of the Hartree-Fock bands of each chain, and the
quartic couplings of the full bond projected onto the Fermi points (C28 eq. 13 convention).

Processes (momentum mismatch modulo 2 pi, 0 = resonant):
  cooper   pair transfer R_A L_A -> R_B L_B between the two right/left pairings, the pairing chosen as
           the one with the smaller mismatch (2 x difference of the pair centres);
  umklapp  R_A R_B <- L_A L_B, mismatch 2 pi (nu - 1) for two crossing bands (nu per cell);
  xback    R_A L_B <- R_B L_A, the inter-pair backscattering (resonant at C28's period-pi reference);
  min_other the smallest mismatch of any other 2 -> 2 process among the four points.
Dimensionless couplings: g_ab = G_ab / (pi (v_a + v_b)) (forward), |p| / (pi (v_A + v_B)) (pair),
|u| / (pi (v_A + v_B)) (umklapp); x, y of C28 eq. (16) for the pair-transfer flow (x < -|y| irrelevant).
usage: resonances.py SPEC [SPEC ...]
Written by Claude (Opus 5.5) for contrib/routeA_n2_chains, 2026-10-28.
"""
import os, sys, json, itertools
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import bonds as B, bands as BD, run_chain as R

d = lambda x: abs((x + np.pi) % (2 * np.pi) - np.pi)


def analyse(spec):
    hspec = spec.replace(':full', ':hf')
    K, meta = R.build(spec)                       # full bond (normalized, fields included)
    Kh, mh = R.build(hspec)
    c, h, U, Hq, H4, r = B.normal_order(Kh)
    T, M = BD.blocks(h)
    pts, _ = BD.fermi_points(T, M, 0.0)
    out = dict(spec=spec, nu_hf=BD.density(T, M, 0.0), n_points=len(pts),
               points=[(round(p['k'], 4), round(p['v'], 4)) for p in pts], dmu=meta['dmu'], scale=meta['scale'])
    if len(pts) != 4:
        return out
    _, _, Uf, _, _, _ = B.normal_order(K)
    Rs = [p for p in pts if p['v'] > 0]; Ls = [p for p in pts if p['v'] < 0]
    best = None
    for perm in ((0, 1), (1, 0)):
        A = (Rs[0], Ls[perm[0]]); Bp = (Rs[1], Ls[perm[1]])
        mm = d(A[0]['k'] + A[1]['k'] - Bp[0]['k'] - Bp[1]['k'])
        if best is None or mm < best[0]:
            best = (mm, A, Bp)
    mm, A, Bp = best
    pair = [A, Bp]
    out['cooper'] = round(mm, 4)
    out['umklapp'] = round(d(Rs[0]['k'] + Rs[1]['k'] - Ls[0]['k'] - Ls[1]['k']), 4)
    out['xback'] = round(d(A[0]['k'] + Bp[1]['k'] - Bp[0]['k'] - A[1]['k']), 4)
    others = []
    for (a, b), (c_, e) in itertools.combinations(list(itertools.combinations(range(4), 2)), 2):
        others.append(d(pts[a]['k'] + pts[b]['k'] - pts[c_]['k'] - pts[e]['k']))
    named = sorted([out['cooper'], out['umklapp'], out['xback']])
    rest = sorted(others)
    for v in named:
        for i, w in enumerate(rest):
            if abs(w - v) < 1e-3:
                rest.pop(i); break
    out['min_other'] = round(rest[0], 4) if rest else None
    out['delta2'] = round(min(d(p['k'] - q['k']) for p, q in itertools.combinations(pts, 2)), 4)
    out['centres'] = [round(float(np.angle(np.exp(1j * (P[0]['k'] + P[1]['k']) / 2))), 4) for P in pair]
    v = [0.5 * (P[0]['v'] - P[1]['v']) for P in pair]
    G = np.array([[BD.amp(Uf, pair[a][0], pair[b][1], pair[a][0], pair[b][1]).real for b in (0, 1)] for a in (0, 1)])
    p = BD.amp(Uf, pair[1][0], pair[1][1], pair[0][0], pair[0][1])
    u = BD.amp(Uf, pair[0][0], pair[1][0], pair[0][1], pair[1][1])
    x = -G[0, 0] / (2 * np.pi * v[0]) - G[1, 1] / (2 * np.pi * v[1]) + (G[0, 1] + G[1, 0]) / (np.pi * (v[0] + v[1]))
    y = abs(p) / np.pi * np.sqrt(1 / (2 * v[0] * v[1]) + 2 / (v[0] + v[1]) ** 2)
    out['v'] = [round(z, 4) for z in v]
    out['g_forward'] = np.round([[G[a, b] / (np.pi * (v[a] + v[b])) for b in (0, 1)] for a in (0, 1)], 4).tolist()
    out['g_pair'] = round(float(abs(p) / (np.pi * (v[0] + v[1]))), 4)
    out['g_umklapp'] = round(float(abs(u) / (np.pi * (v[0] + v[1]))), 4)
    out['x'] = round(float(x), 4); out['y'] = round(float(y), 4)
    return out


if __name__ == '__main__':
    for s in sys.argv[1:]:
        print(json.dumps(analyse(s)), flush=True)

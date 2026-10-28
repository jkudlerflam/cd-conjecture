#!/usr/bin/env python3
"""Tables for the VUMPS (fes_*.json), ring ED (ed_*.json) and open-chain (obc_*.json) runs.

FES: S = mean entropy over the bonds of the unit cell, xi from the transfer matrix (sites).
  c_fit(k) = 6 x slope of S against ln xi over the last k chi values; c_pair between neighbours.
  dS, deb = largest deviation of S and of the bond energy across the unit cell (translation breaking).
ED: rings of L sites (both parity sectors pooled): ground-state energy per site, lowest gap and
  its momentum; L x gap (z = 1 if it converges); the Casimir fit E0/L = e + a/L^2 over the sizes.
OBC: Calabrese-Cardy fit of S(l) = (c/6) ln[(2L/pi) sin(pi l/L)] + s0 over L/4 <= l <= 3L/4,
  after averaging S over a window of 2 pi/(Fermi momentum) sites is not attempted: the raw fit and a
  fit to the running average over 8 consecutive bonds are both reported.
usage: analyze.py fes|ed|obc data/...json
Written by Claude (Opus 5.5) for contrib/routeA_incommensurate, 2026-10-28.
"""
import sys, json
import numpy as np


def fes(path):
    d = json.load(open(path))
    rows = []
    for r in d['results']:
        S = np.array(r['S'])
        rows.append(dict(chi=r['chi'], S=float(S.mean()), dS=float(np.abs(S - S.mean()).max()), xi=r['xi'], e=r['e'],
                         deb=(float(np.ptp(r['ebond'])) if 'ebond' in r else None),
                         dens=np.round(np.mean(r['dens'], axis=1), 4).tolist() if 'dens' in r else None))
    out = dict(path=path, rows=rows)
    if len(rows) >= 3:
        S = np.array([x['S'] for x in rows]); lx = np.log([x['xi'] for x in rows])
        out['c_pair'] = (6 * np.diff(S) / np.diff(lx)).round(3).tolist()
        for k in (3, 4, 5):
            if len(rows) >= k:
                out[f'c_fit{k}'] = round(float(6 * np.polyfit(lx[-k:], S[-k:], 1)[0]), 3)
        out['c_all'] = round(float(6 * np.polyfit(lx, S, 1)[0]), 3)
    return out


def ed(path):
    d = json.load(open(path))
    rows = []
    for r in d['results']:
        L = r['L']
        lev = sorted([(E, k, s['Q']) for s in r['sectors'] for E, k in zip(s['E'], s['k'])])
        E0 = lev[0][0]
        gaps = [(E - E0, k, Q) for (E, k, Q) in lev[1:6]]
        rows.append(dict(L=L, e0=E0 / L, gap=gaps[0][0], Lgap=L * gaps[0][0], k1=round(gaps[0][1] / np.pi, 3),
                         Q0=lev[0][2], gaps=[round(g[0], 6) for g in gaps]))
    out = dict(path=path, rows=rows)
    if len(rows) >= 3:
        L = np.array([x['L'] for x in rows], float); e = np.array([x['e0'] for x in rows])
        A = np.stack([np.ones_like(L), 1 / L ** 2], 1)
        co = np.linalg.lstsq(A, e, rcond=None)[0]
        out['casimir'] = dict(e_inf=float(co[0]), a=float(co[1]), c_v=float(-6 * co[1] / np.pi))
    return out


def obc(path):
    d = json.load(open(path))
    L = d['args']['L']
    r = d['results'][-1]
    S = np.array(r['S'])
    l = np.arange(1, L)
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    m = (l >= L // 4) & (l <= 3 * L // 4)
    c_raw = 6 * np.polyfit(x[m], S[m], 1)[0]
    w = 8
    Sa = np.convolve(S, np.ones(w) / w, mode='valid')
    la = np.arange(len(Sa)) + (w + 1) / 2
    xa = np.log((2 * L / np.pi) * np.sin(np.pi * la / L))
    ma = (la >= L // 4) & (la <= 3 * L // 4)
    c_avg = 6 * np.polyfit(xa[ma], Sa[ma], 1)[0]
    return dict(path=path, L=L, chi=r['chi'], trunc=r['trunc'], c_raw=round(float(c_raw), 3), c_avg=round(float(c_avg), 3))


if __name__ == '__main__':
    mode = sys.argv[1]
    for p in sys.argv[2:]:
        o = dict(fes=fes, ed=ed, obc=obc)[mode](p)
        print(json.dumps(o))

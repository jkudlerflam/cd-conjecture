#!/usr/bin/env python3
"""Finite-entanglement-scaling fits of VUMPS runs: c = 6 dS / d ln(xi) between consecutive bond
dimensions (S averaged over the bonds of the unit cell, xi in logical cells), and a least-squares
slope over the last three points.  usage: fes_fit.py FILE.json [FILE.json ...]
Written by Claude (Opus 5.5) for contrib/routeA_n34_chains, 2026-10-28."""
import sys, json
import numpy as np


def fit(path):
    d = json.load(open(path))
    rs = [r for r in d['results'] if np.isfinite(r['xi']) and 0 < r['xi'] < 1e6]
    chi = [r['chi'] for r in rs]
    S = np.array([np.mean(r['S']) for r in rs]); xi = np.array([r['xi'] for r in rs])
    pairs = [float(6 * (S[i + 1] - S[i]) / np.log(xi[i + 1] / xi[i])) for i in range(len(rs) - 1)]
    last3 = None
    if len(rs) >= 3:
        x = np.log(xi[-3:]); y = S[-3:]
        last3 = float(6 * np.polyfit(x, y, 1)[0])
    m = np.array(chi) >= 16
    c16 = float(6 * np.polyfit(np.log(xi[m]), S[m], 1)[0]) if m.sum() >= 2 else None
    sl_chi = float(np.polyfit(np.log(np.array(chi)[m]), S[m], 1)[0]) if m.sum() >= 2 else None
    return dict(file=path, spec=d['info'].get('spec'), chi=chi, S=S.tolist(), xi=xi.tolist(),
                e=[r['e'] for r in rs], err=[r['err'] for r in rs], c_pairs=pairs, c_last3=last3,
                c_ls_chi16=c16, dS_dlnchi=sl_chi)


if __name__ == '__main__':
    for p in sys.argv[1:]:
        try:
            f = fit(p)
        except Exception as e:
            print(p, 'failed', e); continue
        print(f['spec'])
        print('   chi', f['chi'])
        print('   xi ', np.round(f['xi'], 3).tolist())
        print('   S  ', np.round(f['S'], 4).tolist())
        print('   c (pairs)', np.round(f['c_pairs'], 3).tolist(), ' c(last 3) %s' % (None if f['c_last3'] is None else round(f['c_last3'], 3)),
              ' c(lsq, chi>=16) %s  dS/dln chi %s' % (None if f['c_ls_chi16'] is None else round(f['c_ls_chi16'], 3),
                                                       None if f['dS_dlnchi'] is None else round(f['dS_dlnchi'], 4)))

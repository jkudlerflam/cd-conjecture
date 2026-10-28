#!/usr/bin/env python3
"""Prints the calibration and survey rows quoted in README.md from the data/ files.
Two-size Casimir estimates (L1, L2) are divided by the overshoot of the same two-size fit on
the SU(d) chain at the same sizes (1.053 for d = 4 at L = 8, 12 and 1.117 for d = 5 at L = 5, 10)."""
import glob
import numpy as np
import core
from analyze import load, analyze

R2 = lambda A, hn: A / (2 * np.pi * hn)


def calib():
    rs = load(glob.glob('data/calib.jsonl') + glob.glob('data/big_calib.jsonl'))
    spec = {'XX': dict(Lmin=10), 'Heisenberg': dict(Lmin=10), 'Ising': dict(Lmin=10, nocurrent=True),
            'SU3': dict(Lmod=3), 'SU4': dict(Lmod=4), 'TB': dict(Lmin=8)}
    for (m, p), r in sorted(rs.items()):
        a = analyze(r, **{k: v for k, v in spec[m].items()})
        print('%-11s Ls=%s A_cas=%s c_half=%.3f v(Lmax)=%.4f A_cv=%.4f R_cas=%s R_cv=%.4f' % (
            m, a['Ls'], ('%.4f' % a['A_cas']) if 'A_cas' in a else '-', a['c_half'], a['v_L'][-1],
            a['c_half'] * a['v_L'][-1], ('%.4f' % a['R_cas']) if 'R_cas' in a else '-', R2(a['c_half'] * a['v_L'][-1], a['hnorm'])))


def two_size(path, La, Lb, cal, label, angle=True):
    rs = load([path])
    out = []
    for key, r in rs.items():
        if La not in r or Lb not in r:
            continue
        A = core.fit_cas([La, Lb], [r[La]['E0'], r[Lb]['E0']], quartic=False)[0] / cal
        c = core.c_half(r[La]['S_half'], La, r[Lb]['S_half'], Lb)
        th = (np.arctan2(key[1][1], key[1][0]) % (2 * np.pi)) / np.pi if angle else None
        out.append((th if angle else tuple(np.round(key[1], 3)), r[La]['hnorm'], A, c, R2(A, r[La]['hnorm'])))
    for row in sorted(out, key=lambda x: str(x[0]) if not angle else x[0]):
        print('%-7s %s ||h||=%.4f A=%.4f c_half=%.3f R=%.4f' % (label, ('theta=%.3fpi' % row[0]) if angle else row[0], *row[1:]))


if __name__ == '__main__':
    calib()
    rs = load(['data/blbq_scan.jsonl'])
    for key, r in sorted(rs.items(), key=lambda x: np.arctan2(x[0][1][1], x[0][1][0])):
        a = analyze(r, Lmod=3, Lmin=6)
        print('BLBQ    theta=%.3fpi ||h||=%.4f A_cas=%.4f A_cv=%.4f c_half=%.3f R=%.4f'
              % (np.arctan2(key[1][1], key[1][0]) / np.pi, a['hnorm'], a['A_cas'], a['A_cv'], a['c_half'], a['R_cas']))
    two_size('data/sp4_scan.jsonl', 8, 12, 1.053, 'Sp4')
    two_size('data/sp4pi.jsonl', 8, 12, 1.053, 'Sp4pi', angle=False)
    for f in ('data/so5_scan.jsonl', 'data/so5_neg.jsonl'):
        if glob.glob(f):
            two_size(f, 5, 10, 1.117, 'SO5vec')
    if glob.glob('data/spin32_su4nbhd.jsonl'):
        two_size('data/spin32_su4nbhd.jsonl', 8, 12, 1.053, 'spin3/2', angle=False)

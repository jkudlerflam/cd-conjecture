#!/usr/bin/env python3
"""Markdown tables for README.md from the data files (run analyze.py first or not; this reads data/ directly).

Prints (i) open-chain size series per point, (ii) ring series, (iii) FES fits, (iv) the L=64 cuts,
(v) a power-law extrapolation c_obc(L) = c_inf + a L^-p over the three largest sizes of each liquid point.
Written for contrib/d2_staggered (Claude Opus 5.5 subagent, 2026-10-28).
"""
import os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from scipy.optimize import curve_fit
from analyze import load_obc, fes_fit


def f3(x, nd=3):
    return '-' if x is None else f'{x:.{nd}f}'


def extrap(Ls, cs):
    """c_inf from c = c_inf + a L^-p through the three largest sizes (exact solve when possible)."""
    if len(Ls) < 3:
        return None, None
    L = np.array(Ls[-3:], float); c = np.array(cs[-3:])
    try:
        p, _ = curve_fit(lambda x, ci, a, pw: ci + a * x ** (-pw), L, c, p0=[c[-1], (c[0] - c[-1]) * L[0], 1.0], maxfev=20000)
        if not (0.2 < p[2] < 6):
            raise RuntimeError
        return float(p[0]), float(p[2])
    except Exception:
        # fall back to linear in 1/L through the last two sizes
        ci = c[-1] - (c[-2] - c[-1]) / (L[-1] / L[-2] - 1)
        return float(ci), 1.0


def obc_table(sub='obc', pattern='[st]*'):
    out = ['| model | u | b | L | c_obc (even/odd spread) | c_size | K_F | D2/D1 | x_n (neutral, q=0) | x_c (inside q=1) | max trunc |',
           '| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- | --- | ---: |']
    ext = []
    for g in load_obc(pattern, sub):
        for r in g['rows']:
            out.append('| %s | %s | %s | %d | %.3f (%.3f) | %s | %.2f | %s | %s | %s | %.0e |' % (
                g['name'], g['u'], g['b'], r['L'], r['c_obc'], r['c_obc_spread'], f3(r.get('c_size')), r['K_F'],
                f3(r.get('D2_over_D1')), ', '.join(f'{x:.3f}' for x in r.get('x_n', [])) or '-',
                ', '.join(f'{x:.3f}' for x in r.get('x_c', [])) or '-', r['trunc']))
        Ls = [r['L'] for r in g['rows']]; cs = [r['c_obc'] for r in g['rows']]
        ci, pw = extrap(Ls, cs)
        ext.append((g['name'], g['u'], g['b'], Ls[-1], cs[-1], ci, pw))
    return '\n'.join(out), ext


def fes_table():
    out = ['| model | u | b | chi range | xi at largest chi | c, joint fit of largest four chi | c, consecutive pairs (last three) |',
           '| --- | ---: | ---: | --- | ---: | ---: | --- |']
    for f in sorted(glob.glob(os.path.join(HERE, 'data', 'fes', '*.json'))):
        r = fes_fit(f)
        out.append('| %s | %s | %s | %d to %d | %.1f | %s | %s |' % (
            r['name'], r['u'], r['b'], r['chi'][0], r['chi'][-1], r['xi'][-1], f3(r['c_joint4']),
            ', '.join(f'{p[2]:.3f}' for p in r['c_pairs'][-3:])))
    return '\n'.join(out)


if __name__ == '__main__':
    t, ext = obc_table('obc', '[sti]*')
    print('## open chains\n' + t)
    print('\n## extrapolation c_obc(L) = c_inf + a L^-p, three largest sizes')
    print('| model | u | b | Lmax | c_obc(Lmax) | c_inf | p |\n| --- | ---: | ---: | ---: | ---: | ---: | ---: |')
    for e in ext:
        print('| %s | %s | %s | %d | %.3f | %s | %s |' % (e[0], e[1], e[2], e[3], e[4], f3(e[5]), f3(e[6], 2)))
    t, _ = obc_table('pbc', '*')
    print('\n## rings (c_obc column is the two-endpoint ring fit)\n' + t)
    print('\n## FES\n' + fes_table())
    t, _ = obc_table('obc', 'cut*')
    print('\n## L=64 cuts\n' + t)

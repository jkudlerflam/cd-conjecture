#!/usr/bin/env python3
"""Plot a cut.py output: c_S, K_L, and scaled TRIM gaps against the cut parameter.

usage: cutplot.py data/cuts/X.json [--png out.png]
Scaled gaps L*(E - E0) are shown for the lowest level of the sectors
(0,-), (pi,+), (pi,-) and the first excitation of (0,+).
"""
import sys, json, argparse
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fn'); ap.add_argument('--png', default=None)
    a = ap.parse_args()
    d = json.load(open(a.fn))
    xn = d['x']
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    Ls = sorted({r['L'] for r in d['res']})
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for L in Ls:
        rs = sorted([r for r in d['res'] if r['L'] == L], key=lambda r: r['p'][xn])
        x = np.array([r['p'][xn] for r in rs])
        ax[0, 0].plot(x, [r['cS'] for r in rs], '.-', label=f'L={L}')
        ax[0, 1].plot(x, [r['K'] for r in rs], '.-', label=f'L={L}')
        ax[0, 2].plot(x, [r['mz'] for r in rs], '.-', label=f'L={L}')
        for axx, key, title in ((ax[1, 0], '0,-1', '(0,-)'), (ax[1, 1], f'{L // 2},1', '(pi,+)'),
                                (ax[1, 2], f'{L // 2},-1', '(pi,-)')):
            g = []
            for r in rs:
                lev = r['lev'].get(key, r['lev'].get(key.split(',')[0] + ',0'))
                g.append(L * (lev[0] - r['E0']) if lev else np.nan)
            axx.plot(x, g, '.-', label=f'L={L}'); axx.set_title('L*gap ' + title); axx.set_ylim(0, 8)
    ax[0, 0].set_title('c_S'); ax[0, 0].set_ylim(0, 2); ax[0, 1].set_title('K_L'); ax[0, 1].set_ylim(0, 1.5)
    ax[0, 2].set_title('<Z>')
    for axx in ax.ravel():
        axx.legend(fontsize=7); axx.set_xlabel(xn); axx.grid(alpha=.3)
    fig.suptitle(f"{a.fn} {d['fixed']}")
    fig.tight_layout(); fig.savefig(a.png or a.fn.replace('.json', '.png'), dpi=75)


if __name__ == '__main__':
    main()

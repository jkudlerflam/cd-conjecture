#!/usr/bin/env python3
"""Re-score every logged optimizer point with the drift-penalized two-point estimator and list
the top candidates per d.  Also computes the reference score of the Heisenberg (d = 2) and
SU(3) (d = 3) chains at the same ring sizes.

    python3 rescore.py data/*.jsonl data/typhon/*.jsonl --top 15 --out data/rescore_top.json
"""
import argparse, json, glob
import numpy as np
from core import fit_cas, LOG2


def drift_score(Ls, E0, M, d):
    den = np.pi * M * np.log(d) / LOG2
    two = [fit_cas(Ls[i:i + 2], E0[i:i + 2], quartic=False)[0] / den for i in range(len(Ls) - 1)]
    return min(two) - 0.5 * (max(two) - min(two)), two


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+')
    ap.add_argument('--top', type=int, default=15)
    ap.add_argument('--out', default='data/rescore_top.json')
    a = ap.parse_args()
    rows = []
    for f in a.files:
        for pat in glob.glob(f):
            Ls = None
            for line in open(pat):
                try:
                    r = json.loads(line)
                except Exception:
                    continue
                if 'Ls' in r:
                    Ls = r['Ls']
                info = r.get('info', {})
                if 'E0' not in info or Ls is None or info.get('M', 0) < 1e-9:
                    continue
                d = 2 if ('d2' in pat) else 3 if ('d3' in pat) else int(round(len(r['x']) ** 0.5)) if 'x' in r else 2
                if len(info['E0']) != len(Ls):
                    continue
                sc, two = drift_score(Ls, info['E0'], info['M'], d)
                rows.append({'file': pat, 'd': d, 'start': r.get('start'), 'eval': r.get('eval'), 'final': r.get('final', False),
                             'init': r.get('init'), 'score_new': float(sc), 'two_point': [float(t) for t in two],
                             'ratio_fit': float(fit_cas(Ls, info['E0'])[0] / (np.pi * info['M'] * np.log(d) / LOG2)),
                             'M': info['M'], 'Ls': Ls, 'E0': info['E0'], 'x': r.get('x'), 'h_re': r.get('h_re'), 'h_im': r.get('h_im')})
    out = {}
    for d in sorted({r['d'] for r in rows}):
        rd = sorted([r for r in rows if r['d'] == d], key=lambda r: -r['score_new'])
        out[str(d)] = rd[:a.top]
        print(f'd = {d}: {len(rd)} logged points, {sum(r["final"] for r in rd)} finals')
        for r in rd[:a.top]:
            print('  score %.4f  fit %.4f  two-point %s  M %.3f  %s start %s eval %s %s' % (
                r['score_new'], r['ratio_fit'], [round(t, 3) for t in r['two_point']], r['M'], r['init'], r['start'], r['eval'], r['file'].split('/')[-1]))
    json.dump(out, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()

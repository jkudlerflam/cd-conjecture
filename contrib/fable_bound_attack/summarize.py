#!/usr/bin/env python3
"""Markdown tables of the family scans (data/families_*.jsonl, data/typhon/families_*.jsonl)
and of the optimizer finals, for README.md.  Prints to stdout and writes data/summary.json."""
import json, glob
import numpy as np


def fam_rows(files):
    rows = []
    for f in files:
        for pat in glob.glob(f):
            for line in open(pat):
                try:
                    r = json.loads(line)
                except Exception:
                    continue
                rows.append(r)
    return rows


def fmt(x, n=3):
    return '' if x is None or (isinstance(x, float) and not np.isfinite(x)) else ('%.*f' % (n, x))


def main():
    rows = fam_rows(['data/families_*.jsonl', 'data/typhon/families_*.jsonl'])
    # keep, per (family, par, dens), the run with the most ring sizes (largest L_max)
    best = {}
    for r in rows:
        key = (r['family'], json.dumps(r['par'], sort_keys=True), json.dumps(r.get('dens')))
        if key not in best or (max(r['Ls']), len(r['Ls'])) > (max(best[key]['Ls']), len(best[key]['Ls'])):
            best[key] = r
    out = []
    print('| family | parameters | density | d | rings | M | A (fit, two-size) | R (fit, two-size) | c_half | gap x L |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    for key in sorted(best, key=lambda k: (best[k]['d'], k[0], k[1], k[2])):
        r = best[key]
        gl = r['gapL']
        ok = all(np.isfinite(g) for g in gl)
        line = '| %s | %s | %s | %d | %s | %s | %s, %s | %s, %s | %s | %s |' % (
            r['family'], ', '.join('%s=%s' % (k, v) for k, v in r['par'].items()), ', '.join(fmt(x, 3) for x in (r.get('dens') or [])),
            r['d'], ','.join(map(str, r['Ls'])), fmt(r['M']), fmt(r['A_fit']), fmt(r['A_last2']), fmt(r['ratio']), fmt(r['ratio_last2']),
            fmt(r.get('c_half'), 2), ', '.join(fmt(g, 2) for g in gl) if ok else '')
        print(line)
        out.append({'key': key, 'row': r})
    json.dump(out, open('data/summary.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()

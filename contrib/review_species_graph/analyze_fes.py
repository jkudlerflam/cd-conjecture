#!/usr/bin/env python3
"""Central-charge estimates from the FES runs of fes_effective.py.

c_xi : 6 dS/d ln xi, fitted over the last `k` chi values (and the spread of
       consecutive-pair slopes).
c_chi: from dS/d ln chi = 1/(sqrt(12/c)+1), i.e. c = 12/(1/s - 1)^2 with s the
       fitted slope (Pollmann, Mukerjee, Turner, Moore 2009).
Both are calibrated by the free runs (c = 1 + npay exactly).
usage: analyze_fes.py data/*.json [--k 4]
"""
import sys, json, glob
import numpy as np


def est(path, k=4):
    d = json.load(open(path))
    R = [r for r in d['results'] if r['err'] < 1e-4]
    chi = np.array([r['chi'] for r in R], float)
    S = np.array([np.mean(r['S']) for r in R])
    xi = np.array([r['xi'] for r in R])
    if len(R) < 3:
        return None
    sel = slice(max(0, len(R) - k), len(R))
    sx = np.polyfit(np.log(xi[sel]), S[sel], 1)[0]
    sc = np.polyfit(np.log(chi[sel]), S[sel], 1)[0]
    pairs = 6 * np.diff(S) / np.diff(np.log(xi))
    c_chi = 12 / (1 / sc - 1) ** 2 if 0 < sc < 1 else float('nan')
    return dict(name=d['name'], npay=d['args']['npay'], g=d['args']['g'], seed=d['args']['seed'],
                chis=f"{int(chi[sel][0])}-{int(chi[-1])}", xi_max=float(xi[-1]),
                c_xi=6 * sx, c_xi_pairs=[round(float(x), 3) for x in pairs[-k + 1:]],
                c_chi=c_chi, densities=[round(x, 4) for x in R[-1]['densities']],
                c_free=d['c_free'])


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    k = 4
    rows = [r for p in sorted(sum((glob.glob(a) for a in args), [])) if (r := est(p, k))]
    for r in rows:
        print(f"{r['name']:<22} npay={r['npay']} g={r['g']:<4} chi={r['chis']:<8} xi_max={r['xi_max']:7.1f} "
              f"c_xi={r['c_xi']:.3f} pairs={r['c_xi_pairs']} c_chi={r['c_chi']:.3f} (free {r['c_free']}) dens={r['densities']}")
    json.dump(rows, open('fes_summary.json', 'w'), indent=1)

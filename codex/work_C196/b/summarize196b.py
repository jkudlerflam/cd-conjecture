"""C196b summary (Claude Opus 5.5 agent, 2026-11-10). Usage: python3 summarize196b.py <outdir> [tag ...]
Per chi record of <tag>.jsonl: I1 flag, e, S_bar (mean of the cell-averaged bond entropy of the last two iterations,
read from <tag>.iter.jsonl), xi_even, xi_odd, xi = max, <N_cell>, Var N_cell. Then the local c between consecutive
converged chi and c = 6 x least-squares slope of S_bar against ln xi over the last three converged chi (with xi,
xi_even and xi_odd), and the growth of xi_even, xi_odd over the last chi doubling.
"""
import json
import os
import sys
import glob
import numpy as np


def load(d, tag):
    recs = [json.loads(l) for l in open(os.path.join(d, tag + '.jsonl')) if l.strip()]
    its = [json.loads(l) for l in open(os.path.join(d, tag + '.iter.jsonl')) if l.strip()]
    its = [r for r in its if 'Sbar' in r]
    out = []
    for r in recs:
        last = [x for x in its if x['chi'] == r['chi'] and x['t'] <= r['t']][-2:]
        Sb = float(np.mean([x['Sbar'] for x in last]))
        Nb = float(np.mean([x['Ncell'] for x in last]))
        tm = r['meas'].get('tm', {})
        out.append(dict(chi=r['chi'], I1=r['conv']['I1'], it=r['conv']['iters'], dE1=r['conv']['dE1'],
                        dE2=r['conv']['dE2'], e=r['e_site'], S=Sb, xe=tm.get('xi_even', np.nan),
                        xo=tm.get('xi_odd', np.nan), N=Nb, var=r['varN'], disc=r['disc'], res=r['res_max']))
    return out


def fit(xs, ys):
    A = np.vstack([np.log(xs), np.ones(len(xs))]).T
    return 6 * np.linalg.lstsq(A, np.asarray(ys), rcond=None)[0][0]


def main(d, tags):
    if not tags:
        tags = sorted(os.path.basename(f)[:-6] for f in glob.glob(os.path.join(d, '*.jsonl'))
                      if not f.endswith('.iter.jsonl'))
    for tag in tags:
        rows = load(d, tag)
        print(f'=== {tag}')
        for r in rows:
            r['xi'] = max(r['xe'], r['xo'])
            print(f"chi={r['chi']:4d} I1={str(r['I1']):5s} it={r['it']:3d} dE={r['dE1']:.1e}/{r['dE2']:.1e} "
                  f"e={r['e']:.11f} Sbar={r['S']:.5f} xi_e={r['xe']:.2f} xi_o={r['xo']:.2f} N={r['N']:.4f} "
                  f"varN={r['var']:.4f} disc={r['disc']:.1e} res={r['res']:.1e}")
        conv = [r for r in rows if r['I1']]
        if len(conv) >= 2:
            print('  local c (xi = max):', ' '.join(
                f"{a['chi']}->{b['chi']}:{6 * (b['S'] - a['S']) / np.log(b['xi'] / a['xi']):.3f}"
                for a, b in zip(conv[:-1], conv[1:])))
        if len(conv) >= 3:
            l3 = conv[-3:]
            print(f"  c (last three converged chi {[r['chi'] for r in l3]}): xi=max {fit([r['xi'] for r in l3], [r['S'] for r in l3]):.3f}, "
                  f"xi_even {fit([r['xe'] for r in l3], [r['S'] for r in l3]):.3f}, "
                  f"xi_odd {fit([r['xo'] for r in l3], [r['S'] for r in l3]):.3f}")
            top = conv[-1]
            half = [r for r in conv if r['chi'] * 2 == top['chi']]
            if half:
                h = half[0]
                print(f"  last doubling {h['chi']}->{top['chi']}: xi_even x{top['xe'] / h['xe']:.3f}, "
                      f"xi_odd x{top['xo'] / h['xo']:.3f}")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])

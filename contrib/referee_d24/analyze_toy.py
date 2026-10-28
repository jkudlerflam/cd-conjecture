"""Tables for README.md from results.jsonl (toy_runs.py output).

Writes toy_summary.json.  For each colour number b and ring length L it
reports the grand-canonical ground sector at one fixed chemical potential,
fixed once from the L=12 (b=2) or L=12 (b=3) data, its degeneracy and
momentum, the excess of the best wrong-length (defect) sector, and the
charge-1 and charge-b curvatures.
"""
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
recs = [json.loads(l) for l in open(HERE / 'results.jsonl')]
sec = defaultdict(dict)
for r in recs:
    if r['kind'] == 'sector':
        sec[(r['b'], r['L'])][r['N']] = r

summary = {}
for b in sorted({k[0] for k in sec}):
    if (b, 12) not in sec:
        continue
    E12 = {N: r['E'][0] for N, r in sec[(b, 12)].items()}
    if b == 2 and 4 in E12 and 8 in E12:
        mu = (E12[8] - E12[4]) / 4
    elif 5 in E12 and 7 in E12:
        mu = (E12[7] - E12[5]) / 2
    else:
        continue
    rows = []
    for L in sorted(k[1] for k in sec if k[0] == b):
        d = sec[(b, L)]
        G = {N: r['E'][0] - mu * N for N, r in d.items()}
        Nstar = min(G, key=G.get)
        legal = [N for N in G if N % b == 0]
        wrong = [N for N in G if N % b != 0]
        row = dict(L=L, N_star=Nstar, N_star_legal=(Nstar % b == 0),
                   gs_degeneracy=d[Nstar]['gs_degeneracy'],
                   gs_momentum=d[Nstar]['gs_momenta_2pi_over_L'],
                   gap_in_sector=d[Nstar]['E'][1] - d[Nstar]['E'][0],
                   Ns_computed=sorted(G))
        if legal and wrong:
            row['defect_excess'] = min(G[N] for N in wrong) - min(G[N] for N in legal)
        if all(n in d for n in (Nstar - 1, Nstar, Nstar + 1)):
            row['charge1_curvature'] = (d[Nstar + 1]['E'][0] + d[Nstar - 1]['E'][0]
                                        - 2 * d[Nstar]['E'][0])
        if all(n in d for n in (Nstar - b, Nstar, Nstar + b)):
            row['chargeb_curvature'] = (d[Nstar + b]['E'][0] + d[Nstar - b]['E'][0]
                                        - 2 * d[Nstar]['E'][0])
        rows.append(row)
    summary['b%d' % b] = dict(mu=mu, rows=rows)

dim = [r for r in recs if r['kind'] == 'dimer']
summary['dimer'] = sorted(({k: r[k] for k in ('L', 'J', 'eta', 'alternation', 'zz_AB', 'zz_BA', 'gap')}
                           for r in dim), key=lambda x: (x['L'], x['eta'], x['J']))
(HERE / 'toy_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
for k, v in summary.items():
    print(k)
    if isinstance(v, dict):
        print('  mu =', round(v['mu'], 6))
        for row in v['rows']:
            print('  ', {a: (round(x, 5) if isinstance(x, float) else x) for a, x in row.items()})
    else:
        for row in v:
            print('  ', {a: (round(x, 6) if isinstance(x, float) else x) for a, x in row.items()})

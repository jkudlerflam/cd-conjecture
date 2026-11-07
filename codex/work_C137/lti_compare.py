"""Compare window optima W_n with C135's full LTI values E_LTI(n) (primal objective 'value' in C135 tables.json,
a floating-point SDP value). Reports n^2 (E_LTI(n) - W_n)/M, the loss from restricting X to windows."""
import json, warnings
warnings.filterwarnings('ignore')
from common import C135, HERE
from report import load_survey
T = json.load(open(C135 / 'tables.json'))
rows = []
for name in ['Heisenberg', 'XX', 'XXZ_0.5', 'XXZ_2', 'XXZ_-0.9', 'XXZ_-0.99', 'TFIM_1', 'chiral'] + [f'random_{i:02d}' for i in range(30)]:
    if name not in T: continue
    M = T[name]['M']
    lti = {r['n']: r['value'] for r in T[name]['lti'] if r.get('value') is not None}
    S = {r['n']: r['opt']['value'] for r in load_survey(name)}
    for n in sorted(set(lti) & set(S)):
        rows.append(dict(name=name, n=n, E_LTI=lti[n], W=S[n], loss=n * n * (lti[n] - S[n]) / M))
json.dump(rows, open(HERE / 'results' / 'lti_compare.json', 'w'), indent=1)
import collections
agg = collections.defaultdict(list)
for r in rows: agg[r['name']].append(r)
for k, v in agg.items():
    print(k, ' '.join('n%d:%.4f' % (r['n'], r['loss']) for r in v))
print('max loss', max(r['loss'] for r in rows), 'min loss', min(r['loss'] for r in rows))

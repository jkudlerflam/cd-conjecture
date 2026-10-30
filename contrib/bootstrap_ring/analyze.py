"""Collect out/main/*.json into summary.json and print the tables used in README.md."""
import glob, json, sys
from collections import defaultdict
import numpy as np

from models import KNOWN_A

import re
outs = sys.argv[1:] if len(sys.argv) > 1 else ['out/main']
out = outs[0]
rows = []
for o in outs:
    for f in glob.glob(o + '/*.json'):
        if '/ed_' in f or '/feas' in f or f.endswith('M.json'):
            continue
        r = json.load(open(f))
        r['suffix'] = re.sub(r'.*_s\d+', '', f.rsplit('/', 1)[1][:-5])
        rows.append(r)
ed = {}
for f in glob.glob(out + '/ed_*.json'):
    for r in json.load(open(f)):
        ed[(r['model'], r['L'])] = r['E0']
feas = [json.load(open(f)) for o in outs for f in glob.glob(o + '/feas*.json')]
try:
    Mval = json.load(open(out + '/M.json'))
except FileNotFoundError:
    Mval = {}

EXACT_EINF = {'Heisenberg': 0.25 - np.log(2)}


def fit_inf(Ls, e):
    """Fit e(L) = e_inf + a/L^2 + b/L^4."""
    Ls = np.array(Ls, float); e = np.array(e)
    X = np.vstack([np.ones_like(Ls), Ls ** -2, Ls ** -4]).T
    if len(Ls) < 4:
        X = X[:, :2]
    c, *_ = np.linalg.lstsq(X, e, rcond=None)
    return c


summary = {'ed': {}, 'boot': [], 'feas': feas, 'M': Mval}
models = sorted({m for m, _ in ed})
for m in models:
    Ls = sorted(L for mm, L in ed if mm == m)
    if m.startswith('random'):
        use = [L for L in Ls if L >= 8]
    elif m == 'SU3_swap':
        use = Ls
    else:
        use = [L for L in Ls if L >= 12]
    c = fit_inf(use, [ed[(m, L)] / L for L in use])
    einf = EXACT_EINF.get(m, c[0])
    if m == 'random_real_d2':
        # gapped with a nearly degenerate pair on even rings; odd rings are frustrated
        ev = [L for L in Ls if L % 2 == 0]
        einf = ed[(m, max(ev))] / max(ev); Aed = 0.0
    elif m in EXACT_EINF:
        # A from the exact e_inf at the largest L (Heisenberg has log corrections)
        Lm = max(Ls); Aed = -(6 * Lm ** 2 / np.pi) * (ed[(m, Lm)] / Lm - einf)
    else:
        Aed = -6 * c[1] / np.pi
    summary['ed'][m] = dict(L=Ls, E0=[ed[(m, L)] for L in Ls], e_inf=float(einf), A_fit=float(Aed),
                            A_known=KNOWN_A.get(m))

groups = defaultdict(list)
for r in rows:
    if r.get('value') is None:
        continue
    tag = r['suffix'] or 'full'
    wrap = r['k'] > r['L'] / 2 and tag == '_wrap'
    key = (r['model'], r['k'], r['s'], tag, wrap)
    groups[key].append(r)

for key in sorted(groups):
    m, k, s, tag, wrap = key
    rs = sorted(groups[key], key=lambda r: r['L'])
    einf = summary['ed'].get(m, {}).get('e_inf')
    Aknown = KNOWN_A.get(m)
    tab = []
    for r in rs:
        L = r['L']; lb = r['value']
        E0 = ed.get((m, L))
        item = dict(L=L, bound_per_site=lb, status=r['status'], solver=r['solver'], seconds=r['seconds'])
        if E0 is not None:
            item['E0_per_site'] = E0 / L
            item['gap_total'] = E0 - L * lb
            if Aknown:
                item['gap_over_casimir'] = (E0 - L * lb) / (np.pi * Aknown / (6 * L))
        if einf is not None:
            item['A_cert'] = 6 * L ** 2 / np.pi * (einf - lb)
        tab.append(item)
    fitinfo = None
    Lmin = 14 if s >= 2 else 20
    big = [t for t in tab if t['L'] >= Lmin]
    if not wrap and len(big) >= 3:
        # two-parameter fit l(L)/L = l_inf - a/L^2 with SDP solutions of mixed accuracy
        Ls = np.array([t['L'] for t in big], float); y = np.array([t['bound_per_site'] for t in big])
        X = np.vstack([np.ones_like(Ls), Ls ** -2]).T
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        rms = float(np.sqrt(np.mean((X @ c - y) ** 2)))
        cov = rms ** 2 * np.linalg.inv(X.T @ X) if len(Ls) > 2 else np.zeros((2, 2))
        fitinfo = dict(L_range=[int(Ls.min()), int(Ls.max())], l_inf=float(c[0]),
                       A_relax=float(-6 * c[1] / np.pi), A_relax_err=float(6 / np.pi * np.sqrt(cov[1, 1])),
                       offset=float(einf - c[0]) if einf is not None else None, rms=rms)
    summary['boot'].append(dict(model=m, k=k, s=s, tag=tag, wrap=wrap, rows=tab, fit=fitinfo))

json.dump(summary, open('summary.json', 'w'), indent=1, default=float)

for m, e in summary['ed'].items():
    print('ED', m, 'e_inf=%.6f' % e['e_inf'], 'A_fit=%.4f' % e['A_fit'], 'A_known=', e['A_known'], 'M=', Mval.get(m))
for g in summary['boot']:
    print('\n##', g['model'], 'k=%d s=%d' % (g['k'], g['s']), g['tag'], 'wrap' if g['wrap'] else '', g['fit'])
    for t in g['rows']:
        print('  L=%3d  lb/L=%.7f  E0/L=%s  gap=%s  gap/Casimir=%s  A_cert=%s  %s %s %.0fs' % (
            t['L'], t['bound_per_site'], '%.7f' % t['E0_per_site'] if 'E0_per_site' in t else '-',
            '%.5f' % t['gap_total'] if 'gap_total' in t else '-',
            '%.2f' % t['gap_over_casimir'] if 'gap_over_casimir' in t else '-',
            '%.2f' % t['A_cert'] if 'A_cert' in t else '-', t['status'], t['solver'], t['seconds']))
for f in feas:
    print('FEAS', f)

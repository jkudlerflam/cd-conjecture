"""Markdown tables for README.md from results/n*.json.  Usage: python3 -B make_table.py"""
import json, os, glob
HERE = os.path.dirname(os.path.abspath(__file__))


def e(x, d=1):
    if x is None:
        return '–'
    return ('%.' + str(d) + 'e') % float(x)


rows, rows2 = [], []
for n in range(1, 24):
    f = os.path.join(HERE, 'results', 'n%d.json' % n)
    if not os.path.exists(f):
        continue
    R = json.load(open(f))
    a = R['runs']['C33_kappa1e-3']
    b = R['runs'].get('C33_kappa1/8')
    cs = [float(x['centre_k']) for x in a['bands']]
    du = None if a['dU2'] is None else min(float(a['dU1']), float(a['dU2']))
    du = float(a['dU1']) if du is None else du
    rows.append('| %d | %s | %d | %.2f | %s | %s | %s | %s | %s | %.1f, %.1f | %s | %s |' % (
        n, '–' if a['p'] is None else str(a['p']), a['f'], a['log10_lambda_min_F'],
        e(a['delta2'], 3), e(a['delta4']), e(a['dPT']), e(du), e(max(cs) - min(cs)),
        a['log10_vmin'], a['log10_vmax'],
        e(b['delta4']) if b else '–', e(b['dPT']) if b else '–'))
print('| n | p | f | log10 λ_min(F) | δ₂ | δ₄ | δ_PT | δ_U | centre spread | log10 v_min, v_max | δ₄ at κ = 1/8 | δ_PT at κ = 1/8 |')
print('| --: | --: | --: | --: | --: | --: | --: | --: | --: | :-- | --: | --: |')
print('\n'.join(rows))

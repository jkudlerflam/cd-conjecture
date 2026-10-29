"""Tables from out/*.json: second-order coefficients for many defect patterns, using the saved kernels."""
import json
from pathlib import Path
import numpy as np
from sources import selector, frozen, B_coeff
HERE = Path(__file__).resolve().parent
b = 15; mu = selector(b)


def pats(k):
    g = [b] * k; F = [False]
    P = {}
    for l1, l2 in [(14, 16), (13, 17), (12, 18), (10, 20)]:
        P[f'pair_{l1}_{l2}_adj'] = ([l1, l2] + g[2:], F * k)
        if k >= 4:
            h = k // 2
            P[f'pair_{l1}_{l2}_far'] = ([l1] + g[1:h] + [l2] + g[h + 1:], F * k)
    P['long_30'] = ([30] + g[2:], F * (k - 1))
    if k >= 3: P['long_45'] = ([45] + g[3:], F * (k - 2))
    for l1 in [7, 5, 3, 2]:
        P[f'split_{l1}_{b - l1}'] = ([l1, b - l1] + g[1:], F * (k + 1))
    P['pin'] = (g, [True] + F * (k - 1))
    if k >= 4: P['two_pins_far'] = (g, [True] + F * (k // 2 - 1) + [True] + F * (k - k // 2 - 1))
    return P


def W_of(L, p):
    return sum(l for l in L if l != b) + sum(1 for l, q in zip(L, p) if l == b and q)


res = []
for f in sorted((HERE / 'out').glob('N*_V*.json'), key=lambda f: (json.loads(f.read_text())['density'], json.loads(f.read_text())['N'])):
    d = json.loads(f.read_text()); chi = np.array(d['chi']); N = d['N']; k = N // b
    a = d['diagnostics']['a']
    g = [b] * k; Bg = B_coeff(chi, g, [False] * k, mu); Eg = frozen(g, [False] * k, mu)
    for name, (L, p) in pats(k).items():
        dB = B_coeff(chi, L, p, mu) - Bg; dE = frozen(L, p, mu) - Eg; W = W_of(L, p)
        res.append(dict(N=N, V=d['V'], L=d['L'], density=d['density'], a=a, defect=name, W=W,
                        first_order_over_J=a * dE, delta_B=dB, delta_B_over_W=dB / W,
                        eps_over_Jt=max(0.0, dB) / W, relative_coeff=-dB / (a * dE),
                        sum_rule_error=d['diagnostics']['sum_rule_error']))
(HERE / 'results.json').write_text(json.dumps(res, indent=1) + '\n')
print('| rho | N | L | defect | W | a dE/J | t dB | t dB/W | rel = -t dB/(a dE) |')
print('| --: | --: | --: | :-- | --: | --: | --: | --: | --: |')
for r in res:
    print(f"| {r['density']:.2f} | {r['N']} | {r['L']} | {r['defect']} | {r['W']} | {r['first_order_over_J']:.4e} | "
          f"{r['delta_B']:+.5e} | {r['delta_B_over_W']:+.4e} | {r['relative_coeff']:+.4f} |")

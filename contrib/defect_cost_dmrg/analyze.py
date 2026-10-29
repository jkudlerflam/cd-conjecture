"""Defect costs, residuals and epsilon(J/t) from out/*.json. Writes results.json and prints tables."""
import json, glob, os
import numpy as np
from model import Model

runs = {}
for f in glob.glob('out/*.json'):
    r = json.load(open(f)); runs[r['tag']] = r


def best(r, k=-1):
    return r['hist'][k] if len(r['hist']) >= abs(k) else None


def first_order(r, contact):
    m = Model(r['lengths'], r['pins'], r['V'], r['J'], r['mu'], r['G'])
    return m.first_order(contact)


def W(r, b=15):
    return sum(l for l in r['lengths'] if l != b) + sum(p for l, p in zip(r['lengths'], r['pins']) if l == b)


rows = []
for tag, r in sorted(runs.items()):
    if tag.startswith('free') or tag.startswith('good'):
        continue
    name, N, V = tag.split('_')[:3]
    Jtag = tag.split('_')[3] if len(tag.split('_')) > 3 else None
    g = runs.get(f'good_{N}_{V}_{Jtag}'); fr = runs.get(f'free_{N}_{V}')
    if g is None or fr is None:
        continue
    n = min(len(r['hist']), len(g['hist']))
    if n == 0:
        continue
    contact = fr['hist'][-1]['contact']
    J = r['J']; w = W(r)
    dE = [r['hist'][k]['E'] - g['hist'][k]['E'] for k in range(n)]
    dE1 = first_order(r, contact) - first_order(g, contact)
    dR = [x - dE1 for x in dE]
    chi = r['hist'][n - 1]['chi']
    err = abs(dE[-1] - dE[-2]) if n >= 2 else float('nan')
    rows.append(dict(defect=name, N=r['N'], V=r['V'], L=r['N'] + r['V'], J=J, W=w, chi=chi,
                     cost=dE[-1], cost_chi_change=err, cost_per_JW=dE[-1] / (J * w),
                     first_order_cost=dE1, delta_R=dR[-1], epsilon=max(0.0, -dR[-1]) / (J * w),
                     delta_R_over_J2W=dR[-1] / (J * J * w),
                     free_E_err=fr['hist'][-1]['E'] - fr['E0_exact'], pgmax=max(r['hist'][n - 1]['pgmax'], g['hist'][n - 1]['pgmax']),
                     trunc=max(r['hist'][n - 1]['trunc'], g['hist'][n - 1]['trunc'])))
json.dump(rows, open('results.json', 'w'), indent=1)
print(f"{'defect':8s} {'N':>3s} {'V':>3s} {'J':>7s} {'W':>3s} {'chi':>4s} {'cost':>12s} {'chi_chg':>9s} "
      f"{'cost/JW':>10s} {'dR':>11s} {'eps':>9s} {'dR/J2W':>9s} {'pgmax':>8s}")
for x in sorted(rows, key=lambda x: (x['N'], x['V'], x['defect'], -x['J'])):
    print(f"{x['defect']:8s} {x['N']:3d} {x['V']:3d} {x['J']:7.4f} {x['W']:3d} {x['chi']:4d} {x['cost']:12.5e} "
          f"{x['cost_chi_change']:9.1e} {x['cost_per_JW']:10.3e} {x['delta_R']:11.3e} {x['epsilon']:9.2e} "
          f"{x['delta_R_over_J2W']:9.4f} {x['pgmax']:8.1e}")
print('negative costs:', [(x['defect'], x['N'], x['V'], x['J']) for x in rows if x['cost'] <= 0])

# vacancy-number crossings, N = 30 at J = 0.01: grand potential min_V [E(V) + s V]
cross = []
Vs = sorted({runs[t]['V'] for t in runs if t.startswith('good_N30') and runs[t]['J'] == 0.01})
if len(Vs) >= 3:
    def E(name, V):
        r = runs.get(f'{name}_N30_V{V}_J0.01')
        return r['hist'][-1]['E'] if r and r['hist'] else None
    names = ['good', 'pair', 'long', 'pin']
    tab = {n: {V: E(n, V) for V in Vs} for n in names}
    for V in Vs[:-1]:
        if tab['good'][V] is None or tab['good'][V + 1] is None:
            continue
        s = tab['good'][V] - tab['good'][V + 1]       # good crosses between V and V+1 here
        for n in names[1:]:
            vals = {v: tab[n][v] for v in Vs if tab[n][v] is not None}
            gvals = {v: tab['good'][v] for v in Vs if tab['good'][v] is not None}
            if len(vals) < len(Vs):
                continue
            Gs = min(e + s * v for v, e in vals.items()); Gg = min(e + s * v for v, e in gvals.items())
            cross.append(dict(defect=n, V_cross=(V, V + 1), tilt=s, dG=Gs - Gg,
                              argmin_defect=min(vals, key=lambda v: vals[v] + s * v),
                              canonical_costs={v: vals[v] - gvals[v] for v in Vs}))
    for c in cross:
        print('crossing', c['defect'], c['V_cross'], f"s={c['tilt']:.6f}", f"dG={c['dG']:.4e}",
              'defect argmin V =', c['argmin_defect'], {k: f'{v:.3e}' for k, v in c['canonical_costs'].items()})
    json.dump(cross, open('crossing.json', 'w'), indent=1)

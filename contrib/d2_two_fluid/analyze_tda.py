#!/usr/bin/env python3
"""Summary of the TDA screen.  usage: analyze_tda.py [data/tda_random.json.gz | data/tda_slice.json.gz]
(gzipped lists of screen_tda.py records); writes the survivors to data/<name>_survivors.json."""
import json, gzip, sys
import numpy as np
src = sys.argv[1] if len(sys.argv) > 1 else 'data/tda_random.json.gz'
R = json.load(gzip.open(src, 'rt'))
R = [r for r in R if 'error' not in r]
n = len(R)
def g(r):
    return np.hypot(r['p']['gL'], np.hypot(*r['p']['gR']))
two = [r for r in R if 'two_sea' in r]
print(f"points {n}; two seas {len(two)}")
crit = {
    'pair bound (b.bind>0.05)': lambda r: r['two_sea']['b']['bind'] > 0.05,
    'Qrel=|Q-2k_a|>0.2': lambda r: r['two_sea']['Qrel'] > 0.2,
    'V-mu*>0 (all-up above vacuum)': lambda r: r['p']['V'] - r['mu_star'] > 0,
    'onset density of first sea < 0.25': lambda r: min(r['two_sea']['a']['dens'], r['two_sea']['b']['dens']) < 0.25 and max(r['two_sea']['a']['dens'], r['two_sea']['b']['dens']) < 0.5,
    'all low-order vertex momenta > 0.2': lambda r: r['two_sea']['dmin'] > 0.2,
}
sel = two
for name, f in crit.items():
    sel = [r for r in sel if f(r)]
    print(f"  and {name}: {len(sel)}")
for name, f in crit.items():
    print(f"  {name} alone (of two-sea): {sum(f(r) for r in two)}")
# where: binned by V and g
print('two-sea fraction by V bin:')
for lo in range(-4, 4):
    s = [r for r in R if lo <= r['p']['V'] < lo + 1]
    t = [r for r in s if 'two_sea' in r]
    tb = [r for r in t if r['two_sea']['b']['bind'] > 0.05 and r['two_sea']['Qrel'] > 0.2 and r['p']['V'] - r['mu_star'] > 0]
    print(f"   V in [{lo},{lo+1}): {len(s)} pts, two-sea {len(t)}, bound+Qrel+stable {len(tb)}")
print('by g = |(gL,gR)|:')
for lo in np.arange(0, 4.5, 0.5):
    s = [r for r in R if lo <= g(r) < lo + .5]
    t = [r for r in s if 'two_sea' in r]
    print(f"   g in [{lo:.1f},{lo+.5:.1f}): {len(s)} pts, two-sea {len(t)}")
json.dump([r for r in sel], open(src.replace('.json.gz', '_survivors.json'), 'w'))
for r in sel:
    ts = r['two_sea']; p = r['p']
    print(f"phi={p['phi']:.3f} V={p['V']:+.3f} gL={p['gL']:.3f} gR={np.hypot(*p['gR']):.3f} arg={np.angle(complex(*p['gR'])):+.3f} mu*={r['mu_star']:.3f} dmu={ts['dmu']:.2f} "
          f"a(k={ts['a']['kmin']:+.2f},q={ts['a']['q']:.2f},rho={ts['a']['dens']:.2f}) b(k={ts['b']['kmin']:+.2f},q={ts['b']['q']:.2f},rho={ts['b']['dens']:.3f},bind={ts['b']['bind']:.2f}) dmin={ts['dmin']:.2f} Qrel={ts['Qrel']:.2f}")

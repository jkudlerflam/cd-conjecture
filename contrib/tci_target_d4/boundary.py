#!/usr/bin/env python3
"""Along each eta cut across the left (c = 3 / dimerized) boundary, locate theta* where the
momentum-pi singlet is scale invariant, phi(12) = phi(8) in units 2 pi / L, and interpolate there
c_half(8 -> 12), F/phi and dphi/phi (L = 12 and 8), and phi/current at L = 12.
Tricritical expectation at these sizes (controls/tci_obf.py): F/phi ~ 6.9 and dphi/phi ~ 6.5 at L = 12."""
import glob, sys
import numpy as np
from collections import defaultdict
from analyze import load, levels

files = sys.argv[1:] or glob.glob('data/grid/*.jsonl') + glob.glob('data/fine/*.jsonl') + glob.glob('data/left/*.jsonl')
recs = load(files)


def fp(r):
    L = r['L']
    lv = [x for x in levels(r) if x[0] > 1e-7]
    g = lambda dk, s: [x[0] for x in lv if x[1] == dk and x[2] == s]
    f = lambda l: l[0] if l else np.nan
    return dict(phi=f(g(L // 2, 0)), F=f(g(0, 0)), d=f(g(L // 2 - 1, 0)), cur=f(g(1, 1)), S=r['S'][-1])


cuts = defaultdict(list)
for (t, e), b in recs.items():
    if 8 in b and 12 in b and 0.38 * np.pi < t < 0.5 * np.pi:
        a, c = fp(b[8]), fp(b[12])
        cuts[round(e, 3)].append((t / np.pi, a, c))
print('%6s %6s | %6s | %5s | %5s %5s | %5s %5s | %5s' % ('eta', 'theta', 'p12/p8', 'chalf', 'F/p12', 'd/p12', 'F/p8', 'd/p8', 'p/cur'))
for e in sorted(cuts):
    pts = sorted(cuts[e], key=lambda x: x[0])
    good = [(abs(np.log(c['phi'] / a['phi'])), t, a, c) for t, a, c in pts if np.isfinite(c['phi'] / a['phi'])]
    if not good:
        continue
    # among the scale-invariant points (|ln phi12/phi8| < 0.07) take the softest one, i.e. the one
    # closest to the dimer collapse (smallest phi / current at L = 12)
    good = [g for g in good if g[0] < 0.07]
    if not good:
        continue
    soft = lambda g: g[3]['phi'] / g[3]['cur'] if np.isfinite(g[3]['cur']) else g[3]['phi']
    dev, t, a, c = min(good, key=soft)
    print('%6.2f %6.3f | %6.3f | %5.2f | %5.2f %5.2f | %5.2f %5.2f | %5.2f' % (
        e, t, c['phi'] / a['phi'], 3 * (c['S'] - a['S']) / np.log(1.5), c['F'] / c['phi'], c['d'] / c['phi'],
        a['F'] / a['phi'], a['d'] / a['phi'], c['phi'] / c['cur']))


# ---------------------------------------------------------------- table and figure
rows = []
for e in sorted(cuts):
    pts = sorted(cuts[e], key=lambda x: x[0])
    good = [(abs(np.log(c['phi'] / a['phi'])), t, a, c) for t, a, c in pts if np.isfinite(c['phi'] / a['phi'])]
    good = [g for g in good if g[0] < 0.07]
    if not good:
        continue
    soft = lambda g: g[3]['phi'] / g[3]['cur'] if np.isfinite(g[3]['cur']) else g[3]['phi']
    dev, t, a, c = min(good, key=soft)
    rows.append(dict(eta=e, theta=t, ratio=c['phi'] / a['phi'], chalf=3 * (c['S'] - a['S']) / np.log(1.5),
                     F12=c['F'] / c['phi'], d12=c['d'] / c['phi'], F8=a['F'] / a['phi'], d8=a['d'] / a['phi'],
                     soft=c['phi'] / c['cur']))
import json
json.dump(rows, open('data/boundary_table.json', 'w'), indent=0)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
E = [r['eta'] for r in rows]
ax[0].plot(E, [r['F12'] for r in rows], 'o-', label='F/phi, L = 12')
ax[0].plot(E, [r['d12'] for r in rows], 's-', label='desc/phi, L = 12')
ax[0].plot(E, [r['F8'] for r in rows], 'o:', c='C0', alpha=0.5, label='F/phi, L = 8')
ax[0].plot(E, [r['d8'] for r in rows], 's:', c='C1', alpha=0.5, label='desc/phi, L = 8')
ax[0].axhspan(6.3, 8.1, color='C3', alpha=0.15, label='tricritical control, L = 8 to 16')
ax[0].axhline(2, color='k', lw=0.8, ls='--', label='free Majorana (descendant)')
ax[0].set_xlabel('eta'); ax[0].set_ylabel('ratio'); ax[0].set_ylim(0, 9); ax[0].legend(fontsize=7)
ax[0].set_title('softest scale-invariant pi singlet on each eta cut')
ax[1].plot(E, [r['chalf'] for r in rows], 'o-')
ax[1].axhspan(2.9, 3.1, color='C0', alpha=0.15, label='c = 3 calibration band')
ax[1].axhline(3.2, color='C3', ls='--', label='16/5')
ax[1].set_xlabel('eta'); ax[1].set_ylabel('c_half(8 -> 12)'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('boundary.png', dpi=110)

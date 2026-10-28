#!/usr/bin/env python3
"""Phase map of the spin-3/2 chart from chart.py output (L = 8, 12).
Classes: FM (ground spin 3L/2), PP (0 < ground spin < 3L/2), D (dimer doublet: the k0+pi singlet
collapses, pi0(12) < 0.05 and pi0(12) < 0.5 pi0(8) in units 2pi/L), C3 (k0 = pi at L = 12, spin 0,
not D), C1 (k0 = 0 at L = 12, lowest excitation a spin-1 level at k0+pi), X (anything else).
Writes map_table.json and map.png."""
import json, glob, sys
import numpy as np
from analyze import load, summary

files = sys.argv[1:] or glob.glob('data/grid/*.jsonl') + glob.glob('data/fine/*.jsonl') + glob.glob('data/left/*.jsonl')
recs = load(files)
rows = []
for key, b in recs.items():
    if 8 not in b or 12 not in b:
        continue
    s8, s12 = summary(b[8]), summary(b[12])
    ch = 3 * (s12['S_half'] - s8['S_half']) / np.log(1.5)
    S = s12['Sgs']
    if S == 18:
        cl = 'FM'
    elif S > 0:
        cl = 'PP'
    elif s12['pi0'] < 0.05 and s12['pi0'] < 0.5 * s8['pi0']:
        cl = 'D'
    elif s12['k0'] == 6:
        cl = 'C3'
    elif s12['k0'] == 0 and s12['gap_dk'] == 6 and s12['gap_S'] == 1:
        cl = 'C1'
    else:
        cl = 'X'
    z = -np.log((s12['gap'] / 12) / (s8['gap'] / 8)) / np.log(1.5) if s12['gap'] > 0 and s8['gap'] > 0 else np.nan
    rows.append(dict(theta=key[0] / np.pi, eta=key[1], cls=cl, chalf=ch, z=z, k0_12=s12['k0'], S12=S,
                     pi0_8=s8['pi0'], pi0_12=s12['pi0'], cur_8=s8['cur'], cur_12=s12['cur'], k00_12=s12['k00'],
                     gap_12=s12['gap'], gap_dk=s12['gap_dk'], gap_S=s12['gap_S']))
json.dump(rows, open('data/map_table.json', 'w'), indent=0)
from collections import Counter
print(Counter(r['cls'] for r in rows))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
mk = dict(FM=('s', '#888888'), PP=('D', '#bbbbbb'), D=('x', '#d62728'), C3=('o', '#1f77b4'), C1=('^', '#2ca02c'), X=('.', '#9467bd'))
for cl, (m, col) in mk.items():
    rr = [r for r in rows if r['cls'] == cl]
    ax[0].scatter([r['theta'] for r in rr], [r['eta'] for r in rr], marker=m, c=col, s=18, label=cl)
ax[0].set_xlabel('theta / pi'); ax[0].set_ylabel('eta'); ax[0].legend(fontsize=8, ncol=3)
ax[0].set_title('phase classes (L = 8, 12)')
rr = [r for r in rows if r['cls'] in ('C3', 'X', 'C1', 'D')]
sc = ax[1].scatter([r['theta'] for r in rr], [r['eta'] for r in rr], c=[min(max(r['chalf'], 1.0), 3.6) for r in rr], cmap='viridis', s=22)
plt.colorbar(sc, ax=ax[1], label='c_half(8 -> 12)')
ax[1].set_xlabel('theta / pi'); ax[1].set_title('c_half, non-ferromagnetic points')
plt.tight_layout(); plt.savefig('map.png', dpi=110)

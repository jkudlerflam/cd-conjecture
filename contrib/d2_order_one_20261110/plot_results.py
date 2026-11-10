#!/usr/bin/env python3
"""Static research figure. Run analyze.py after all selected sectors finish."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/d2-order-one-matplotlib')
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
labels={
 'bridge_r1.5_u1.5_k3':'Ordered flag (R=2)',
 'bridge_r1.5_u2_k1':'Initial flag (R=1.75)',
 'bridge_r2_u1.5_k4':'Simple seed (R=2.67)',
 'bridge_r3_u1.5_k4':'Nearer bridge (R=9.33)',
 'bridge_r3_u1_k4':'Nearer bridge (R=14)',
 'bridge_r4_u1_k4':'Bridge (R=28)',
 'attractive_Jz-0.5':'Attractive rational (R=2)',
 'reference_peotta_exploratory':'Attractive reference (R=7.77)',
}


def main():
    rows=json.loads((HERE/'dmrg_summary.json').read_text())
    rows={r['name']:r for r in rows}
    fig,axs=plt.subplots(1,2,figsize=(12,5.8),layout='constrained')
    ax=axs[0]
    colors=plt.get_cmap('tab10').colors
    for idx,(name,label) in enumerate(labels.items()):
        if name not in rows:continue
        for p in rows[name]['pairs']:
            offset=0 if p['sizes'][0]==24 else .2
            y=idx+offset
            x=[p['c_central3'],p['c_fractional']]
            ax.plot(x,[y,y],color=colors[idx],lw=2)
            marker='o' if p['sizes'][0]==24 else 's'
            ax.scatter(x,[y,y],marker=marker,color=colors[idx],s=28,zorder=4)
    ax.axvline(1,color='.5',ls='--',lw=1);ax.axvline(2,color='.5',ls=':',lw=1)
    ax.set_yticks(range(len(labels)),labels.values(),fontsize=9)
    ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    ax.set_xlabel('Entropy size slope; segment joins the two windows')
    ax.set_title('Circles: 24 → 48 sites   Squares: 48 → 96',fontsize=11)
    ax=axs[1]
    for name,color in [('bridge_r2_u1.5_k4',colors[2]),('bridge_r4_u1_k4',colors[5])]:
        if name not in rows:continue
        ss=[s for s in rows[name]['states'] if s['parity_complete']]
        if not ss:continue
        s=max(ss,key=lambda s:s['L'])
        p=HERE/'dmrg'/f"{name}_L{s['L']}_P{s['parity']}_s{s['seed']}.json"
        r=json.loads(p.read_text());n=np.asarray(r['measurements']['dens'])
        ax.plot((np.arange(len(n))+.5)/len(n),n,'o-',ms=2,lw=.8,color=color,label=labels[name]+f", L={len(n)}")
    ax.set_xlabel('Position / chain length');ax.set_ylabel('Particle density (1 − ⟨Z⟩)/2')
    ax.set_title('Bulk oscillations remain visible',fontsize=11)
    ax.legend(fontsize=8,loc='best');ax.grid(alpha=.15)
    fig.suptitle('Comparable-coupling search: finite-size diagnostics, not CFT determinations',fontsize=13)
    fig.savefig(HERE/'diagnostics.png',dpi=180)
    fig.savefig(HERE/'diagnostics.svg')
    svg=HERE/'diagnostics.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')


if __name__=='__main__':main()

#!/usr/bin/env python3
"""Plot completed, parity-selected data; no active MPS is read."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/d2-order-one-matplotlib')
import json,hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
LABELS={'attr_ref':'Reference R=7.77','pref_ref':'Preferred reference R=2.66','pref_rat':'Rational (14,-5,11,-12), R=2.8','pref_zrat':'Rational (7,-3,6,-9), R=3'}
COLORS={'attr_ref':'#286ca8','pref_ref':'#b54c39','pref_rat':'#8054a0','pref_zrat':'#b37a1a'}

def main():
    s=json.loads((HERE/'summary.json').read_text())
    d=json.loads((HERE.parent/'diagnostics/attractive_diagnostics.json').read_text())
    fig,ax=plt.subplots(1,3,figsize=(13,4.3),layout='constrained')
    for n,label in LABELS.items():
        color=COLORS[n];pairs=[p for p in s['pairs'] if p['name']==n]
        if pairs:
            x=[p['sizes'][1] for p in pairs]
            ax[0].plot(x,[p['c_fractional'] for p in pairs],'o-',color=color,label=label)
            ax[0].plot(x,[p['c_central3'] for p in pairs],'s:',color=color,ms=4)
        selected=s['selected'].get(n,{})
        ll=sorted(map(int,selected))
        if ll:ax[1].plot(ll,[selected[str(L)]['parity_gap'] for L in ll],'o-',color=color)
        order=[]
        for L in ll:
            P=selected[str(L)]['parity']
            rr=[r for r in d['rows'] if r['name']==n and r['L']==L and r['parity']==P]
            if not rr:continue
            r=rr[-1];src=ROOT/r['source_record']
            if hashlib.sha256(src.read_bytes()).hexdigest()!=r['record_sha256']:continue
            key='Z_zero' if n=='pref_zrat' else 'X_pi'
            order.append((L,r['order'][key]['order_square']))
        if order:ax[2].plot(*zip(*order),'o-',color=color)
    ax[0].axhline(1,c='.5',lw=1,ls='--');ax[0].set_ylabel('Entropy size slope')
    ax[0].set_title('Solid: fractional window; dotted: central 3')
    ax[0].set_xlabel('Upper length in the doubled-size pair')
    ax[0].set_ylim(.85,1.55)
    ax[1].set_yscale('log');ax[1].set_ylabel('Parity energy splitting')
    ax[1].set_title('Parity splitting');ax[1].set_xlabel('Chain length')
    ax[2].set_ylabel(r'Central-window order $\langle Q^2\rangle/W^2$')
    ax[2].set_title('Staggered X; uniform Z for the R=3 reserve')
    ax[2].set_xlabel('Chain length')
    for a in ax:
        a.set_xticks([48,96,192]);a.grid(alpha=.18);a.tick_params(labelsize=9)
        a.title.set_fontsize(10)
    fig.legend(handles=[Line2D([0],[0],color=COLORS[n],marker='o',label=label) for n,label in LABELS.items() if s['selected'].get(n)],loc='outside lower center',ncol=2,frameon=False,fontsize=9)
    fig.suptitle('Attractive-chain diagnostics: completed parity-selected states only',fontsize=13)
    fig.savefig(HERE/'diagnostics.png',dpi=180);fig.savefig(HERE/'diagnostics.svg')
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Static scientific summary of every complete full-DMRG field point."""
import json
import os
from fractions import Fraction
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/d2-bridge-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME','/private/tmp/d2-bridge-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from summarize import HERE, summarize

out=summarize()
rows=[r for r in out['rows'] if r['L']==24 and not r['name'].startswith('control')]
families=sorted({r['parameters']['family'] for r in rows})
fig,axes=plt.subplots(3,len(families),figsize=(13,9),sharex='col')
for j,family in enumerate(families):
    for q,marker in [('0','o'),('1/16','s'),('1/8','^')]:
        rr=sorted([r for r in rows if r['parameters']['family']==family and r['parameters']['q']==q],key=lambda r:Fraction(r['parameters']['hz']))
        h=[float(Fraction(r['parameters']['hz'])) for r in rr]
        for i,key in enumerate(['density','S_fractional','density_bulk_range']):
            axes[i,j].plot(h,[r[key] for r in rr],marker=marker,label=f'q={q}')
            axes[i,j].grid(alpha=.25)
    axes[0,j].axhline(.275,ls=':',color='gray',lw=1)
    parameters=next(r['parameters'] for r in rows if r['parameters']['family']==family)
    axes[0,j].set_title(f"r = {parameters['r']}, D = {parameters['D']}")
    axes[0,j].legend(frameon=False)
    axes[2,j].set_xlabel('physical field h')
axes[0,0].set_ylabel('Central-half density')
axes[1,0].set_ylabel('Fractional-window entropy')
axes[2,0].set_ylabel('Central-half density range')
fig.suptitle('Full microscopic L=24 DMRG, lower of both exact parity sectors\nLines guide the eye; raw density includes virtual dressing')
fig.tight_layout()
fig.savefig(HERE/'field_scan.png',dpi=180)
plt.close(fig)

pairs=out['pairs']
if pairs:
    fig,ax=plt.subplots(figsize=(10,max(3,.5*len(pairs))))
    for i,p in enumerate(pairs):
        ax.plot([p['c_fractional'],p['c_central3']],[i,i],'-',color='gray',lw=1)
        ax.plot(p['c_fractional'],i,'o',color='C0',label='fractional window' if i==0 else None)
        ax.plot(p['c_central3'],i,'x',color='C1',label='central three' if i==0 else None)
    ax.set_yticks(range(len(pairs)),[p['name']+' '+str(p['sizes']) for p in pairs],fontsize=8)
    ax.axvline(1,ls=':',color='gray')
    ax.set_xlabel('Two-size entropy slope 6 ΔS / log(L₂/L₁)')
    ax.set_title('Screening observables; finite-size slopes are not CFT identifications')
    ax.grid(axis='x',alpha=.25);ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(HERE/'size_slopes.png',dpi=180);plt.close(fig)

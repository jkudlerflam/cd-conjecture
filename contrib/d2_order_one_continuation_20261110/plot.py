#!/usr/bin/env python3
"""Scientific summary of completed bridge measurements; no active sweeps."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/d2-order-one-continuation-mpl')
os.environ.setdefault('XDG_CACHE_HOME','/private/tmp/d2-order-one-continuation-cache')
from pathlib import Path
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

HERE=Path(__file__).resolve().parent
NAME='b_r5_2_d1_10_q0_hm7_8'

def main():
    summary=json.loads((HERE/'size_comparisons.json').read_text())[NAME]
    gs=[r for r in summary['ground_branches'] if r['both_parities_complete'] and r['convergence']['sweep_converged']]
    assert len(gs)>=3
    ls=np.array([r['L'] for r in gs]);sf=np.array([r['S_fractional'] for r in gs])
    sc=np.array([r['S_central3'] for r in gs]);reference=next(r for r in gs if r['L']==48)['S_fractional']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'axes.titleweight':'bold','axes.labelcolor':'#253449',
        'text.color':'#253449','xtick.color':'#526174','ytick.color':'#526174'})
    fig,axes=plt.subplots(1,3,figsize=(13.8,4.7),layout='constrained')
    blue='#1766a4';purple='#7752a3';orange='#c47626'
    ax=axes[0]
    guide=np.geomspace(24,max(192,int(ls.max())),200)
    for c,style,color in [(1,'--','#98a3ad'),(1.5,':','#697783'),(2,'-.','#c0c6cb')]:
        ax.plot(guide,c*np.log(guide/48)/6,style,color=color,lw=1.2,label=f'c = {c:g} guide')
    ax.plot(ls,sf-reference,'o-',color=blue,lw=2,label='Measured entropy')
    ax.set(xscale='log',xlabel='Chain length L',ylabel='Fractional entropy − S(48)',title='Fixed Hamiltonian, growing chains')
    ax.set_xticks([24,48,64,96,192],[24,48,64,96,192]);ax.grid(alpha=.15)
    ax.legend(fontsize=8,loc='upper left',frameon=False)
    ax=axes[1]
    x=np.sqrt(ls[:-1]*ls[1:]);den=np.log(ls[1:]/ls[:-1])
    ax.axhline(1,color='#98a3ad',ls='--',lw=1.2,label='c = 1 reference')
    ax.plot(x,6*np.diff(sf)/den,'o-',color=blue,lw=2,label='Fractional window')
    ax.plot(x,6*np.diff(sc)/den,'s-',color=purple,lw=1.7,label='Central three cuts')
    ax.set(xscale='log',xlabel='Geometric mean of adjacent lengths',ylabel='Effective entropy coefficient',title='Adjacent-size consistency')
    ax.set_xticks([32,48,64,96,144],[32,48,64,96,144]);ax.set_ylim(.9,1.85);ax.grid(alpha=.15)
    ax.legend(fontsize=8,loc='upper right',frameon=False)
    ax=axes[2]
    records=json.loads((HERE/'diagnostics/bridge_diagnostics.json').read_text())['rows']
    dr=[]
    for g in gs:
        matches=[r for r in records if r['L']==g['L'] and r['parity']==g['parity']]
        if matches:dr.append(max(matches,key=lambda r:r['chi']))
    xx=np.array([r['L'] for r in dr])
    series=[('Transverse window m²',blue,[r['order']['X_pi']['connected_order_square'] for r in dr]),
            ('Far transverse RMS',purple,[r['envelope']['X']['far']['rms'] for r in dr]),
            ('Density modulation RMS',orange,[r['density']['rms_detrended'] for r in dr])]
    for label,color,values in series:
        y=np.array(values);ax.plot(xx,y/y[0],'o-',color=color,lw=1.8,label=label)
    ax.set(xscale='log',yscale='log',xlabel='Chain length L',ylabel='Value / value at smallest length',title='Order and density diagnostics')
    ax.set_xticks([24,48,96,192],[24,48,96,192]);ax.grid(alpha=.15,which='both')
    ax.legend(fontsize=8,loc='lower left',frameon=False)
    ymin=min(min(np.array(values)/values[0]) for _,_,values in series)
    yticks=[t for t in (.01,.02,.05,.1,.2,.5,1.) if t>=ymin*.85]
    ax.set_yticks(yticks,[f'{t:g}' for t in yticks]);ax.yaxis.set_minor_locator(NullLocator())
    for ax in axes:ax.xaxis.set_minor_locator(NullLocator())
    fig.suptitle('Nearest-neighbour qubit lead: Jx = 1/2, Jy = −1/2, Jz = 5/8, D = 1/10, hz = −7/8',fontsize=13,fontweight='bold')
    fig.supxlabel('Completed states only. Finite-size evidence is not a CFT identification; small periodic rings remain close to c = 1.',fontsize=9,color='#566675')
    fig.savefig(HERE/'candidate_diagnostics.png',dpi=160,bbox_inches='tight',pad_inches=.16)
    fig.savefig(HERE/'candidate_diagnostics.svg',bbox_inches='tight',pad_inches=.16)
    plt.close(fig)

if __name__=='__main__':main()

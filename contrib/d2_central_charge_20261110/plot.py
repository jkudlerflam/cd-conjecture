"""Scientific comparison of fixed-Hamiltonian finite-size and bulk tests."""
import os
os.environ['MPLCONFIGDIR']='/private/tmp/d2-central-charge-mplconfig'
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator,ScalarFormatter
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'summary.json').read_text())
old=json.loads((HERE/'inputs/size_comparisons.json').read_text())['b_r5_2_d1_10_q0_hm7_8']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold'})
fig,grid=plt.subplots(2,2,figsize=(12,8.1),constrained_layout=True)
axs=grid.ravel()
ax=axs[0]
rr=[r for r in old['pairs'] if r['sizes'][1]==2*r['sizes'][0]]
ax.plot([np.sqrt(np.prod(r['sizes'])) for r in rr],[r['c_fractional'] for r in rr],'o-',label='Open chain',color='#2468a0')
rr=[r for r in d['periodic_pairs'] if r['model']=='candidate' and min(r['L'])>=16]
if rr:ax.plot([np.sqrt(np.prod(r['L'])) for r in rr],[r['c']['1']['fractional'] for r in rr],'s-',label='Periodic ring',color='#c66a22')
for c in [1.,1.5]:ax.axhline(c,color='0.6',ls='--',lw=1);ax.text(17,c+.015,f'{c:g} guide',color='0.45',fontsize=8)
ax.set(xscale='log',xlabel='Geometric mean of chain lengths',ylabel='Effective entropy coefficient',title='Finite-size entropy')
ax.set_xticks([16,32,64,128,192]);ax.xaxis.set_major_formatter(ScalarFormatter());ax.xaxis.set_minor_locator(NullLocator());ax.legend(loc='upper right')
ax=axs[1]
families=[r for r in d['infinite'] if r['cell']==2 and not r['one_site'] and ((r['model']=='candidate' and r['gauge']) or (r['model']=='xx' and r['seed']==23))]
for r in families:
    pts=r['points'];label='Candidate, real frame' if r['model']=='candidate' else 'XX control'
    color='#2468a0' if r['model']=='candidate' else '#6b8e50'
    ax.plot([p['xi'] for p in pts],[p['S'] for p in pts],'o-',label=label,color=color)
    annotate={16,24,48,96,max(p['chi'] for p in pts)}
    for p in pts:
        if p['chi'] in annotate:ax.annotate(str(p['chi']),(p['xi'],p['S']),xytext=(3,5),textcoords='offset points',fontsize=8,color=color)
ax.set(xscale='log',xlabel='Transfer correlation length ξ (sites)',ylabel='Mean half-chain entropy',title='Infinite-chain entanglement')
limits=ax.get_xlim();ax.set_xticks([25,50,100,200,400,800]);ax.set_xlim(limits)
ax.xaxis.set_major_formatter(ScalarFormatter());ax.xaxis.set_minor_locator(NullLocator())
ax.legend(fontsize=8);ax.grid(alpha=.15)
ax=axs[2]
for r in families:
    if r['model']=='candidate':
        rr=[x for x in r['pairs'] if min(x['chis'])>=48]
        for alpha,color in [('1','#2468a0'),('2','#c66a22'),('3','#8454a1')]:
            ax.plot([np.sqrt(np.prod(x['chis'])) for x in rr],[x['c_renyi'][alpha] for x in rr],'o-',color=color,label=f'Candidate, Rényi order {alpha}')
    else:
        rr=[x for x in r['pairs'] if min(x['chis'])>=24]
        ax.plot([np.sqrt(np.prod(x['chis'])) for x in rr],[x['c'] for x in rr],'s--',color='#6b8e50',label='XX control, order 1')
for c in [.5,1.]:ax.axhline(c,color='0.65',ls=':',lw=1)
for model,label,color,marker in [('ising_critical','Critical Ising','#9a6c2b','^'),('ising_gapped','Gapped Ising','#777777','v')]:
    matches=[r for r in d['infinite'] if r['model']==model and r['cell']==2 and r['seed']==23]
    for r in matches:
        rr=r['pairs']
        if rr:
            update='one-site' if r['one_site'] else 'two-site'
            ax.plot([np.sqrt(np.prod(x['chis'])) for x in rr],[x['c'] for x in rr],marker+'--',color=color,label=f'{label}, {update}')
ax.set(xlabel='Geometric mean of bond dimensions',ylabel='Effective coefficient from successive states',title='Successive bulk entropy coefficients',ylim=(-.03,1.16))
ax.grid(alpha=.15);ax.legend(fontsize=8,loc='upper right')
ax=axs[3]
for r in families:
    if r['model']!='candidate':continue
    pts=r['points'];order=[np.sqrt(np.mean(np.array(p['X'])**2+np.array(p['Y'])**2)) for p in pts]
    ax.plot([p['chi'] for p in pts],order,'o-',color='#8454a1',label='Transverse one-point order')
ax.set(xlabel='Bond dimension χ',ylabel=r'$\sqrt{\overline{\langle X\rangle^2+\langle Y\rangle^2}}$',title='Finite-χ branch change',ylim=(-.01,.6))
ax.grid(alpha=.15);ax.legend(fontsize=8)
fig.suptitle('Fixed qubit Hamiltonian: testing the central-charge interpretation',fontsize=15,fontweight='bold')
fig.text(.5,-.035,'Accepted states only. These finite-scale coefficients do not by themselves identify a thermodynamic central charge.',ha='center',fontsize=10,color='0.35')
for ext in ['png','svg']:
    path=HERE/f'central_charge_diagnostics.{ext}';fig.savefig(path,dpi=160,bbox_inches='tight',pad_inches=.16)
    if ext=='svg':path.write_text('\n'.join(x.rstrip() for x in path.read_text().splitlines())+'\n')

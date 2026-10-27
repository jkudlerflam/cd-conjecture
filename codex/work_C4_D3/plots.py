#!/usr/bin/env python3
import os
from pathlib import Path
Path('mplconfig_d3').mkdir(exist_ok=True)
os.environ['MPLCONFIGDIR']=str(Path('mplconfig_d3').resolve())
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from analyze import load,cf
rows=load();plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':140})
fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
ax=axes[0,0]
for L in [6,8,10]:
 rs=sorted([r for r in rows if r['source']=='scan_sp4.jsonl' and r['L']==L],key=lambda r:r['label'])
 ax.plot([int(r['label'][-3:])/48 for r in rs],[cf(r) for r in rs],'.-',lw=.8,label=f'L={L}')
ax.axhline(3,color='black',ls='--',lw=1);ax.axvspan(1.25,2,color='gray',alpha=.15);ax.set(xlim=(0,2),ylim=(0,10),xlabel=r'$\theta/\pi$',ylabel=r'$c_{\rm ent}$',title='Sp(4): raw entropy slopes');ax.legend(ncol=3,fontsize=8)
ax=axes[0,1]
for t in [.015,.06544984695,.1,.16,.1963495408,.2617993878]:
 lab=f'sp4_angle{t:.11f}';rs=sorted([r for r in rows if r['label']==lab and r['L']>=8],key=lambda r:r['L'])
 ax.plot([r['L'] for r in rs],[r['v_current_assumed'] for r in rs],'o-',label=f'{t:.3f}')
ax.set(xlabel='L',ylabel='Adjacent-gap velocity',title='Sp(4): slow-mode candidates');ax.legend(title=r'$\theta$',ncol=2,fontsize=8)
ax=axes[0,2];rs=[r for r in rows if r['source']=='scan_su2.jsonl' and r['L']==10]
i=np.array([int(r['label'][-3:]) for r in rs]);z=1-2*(i+.5)/192;phi=(i*np.pi*(3-np.sqrt(5)))%(2*np.pi)
sc=ax.scatter(phi/np.pi,z,c=[cf(r) for r in rs],vmin=0,vmax=6,cmap='viridis',s=35);fig.colorbar(sc,ax=ax,label=r'$c_{\rm ent}(L=10)$');ax.set(xlabel=r'$\phi/\pi$',ylabel=r'$\cos\vartheta$',title='Spin 3/2: projector sphere')
ax=axes[1,0];vals=[-2,-1,-.75,-.5,-.25,0,.125,.25,.375,.5,.75,1,2];arr=np.full((len(vals),len(vals)),np.nan)
for r in rows:
 if r['source']=='scan_orbital_dense.jsonl' and r['L']==8 and r['coeff'][0]==1:
  x=r['coeff'][2];y=r['coeff'][1];i=vals.index(x);j=vals.index(y);arr[i,j]=arr[j,i]=cf(r)
im=ax.imshow(arr,origin='lower',vmin=0,vmax=6,cmap='viridis',aspect='auto');ax.set_xticks(range(len(vals)),[str(v) for v in vals],rotation=90,fontsize=7);ax.set_yticks(range(len(vals)),[str(v) for v in vals],fontsize=7);ax.set(xlabel='y',ylabel='x',title='Spin-orbital + product, L=8');fig.colorbar(im,ax=ax,label=r'$c_{\rm ent}$')
ax=axes[1,1]
for L in [6,9,12]:
 rs=sorted([r for r in rows if r['family']=='spin1' and r['label'].startswith('bb_t') and r['L']==L and r['source'] in ['scan_spin1.jsonl','corrected_spin1_L9.jsonl']],key=lambda r:r['label'])
 ax.plot([int(r['label'][-3:])/48 for r in rs],[cf(r) for r in rs],'.-',lw=.8,label=f'L={L}')
ax.axhline(2,color='black',ls='--',lw=1);ax.axvspan(.5,1.25,color='gray',alpha=.15);ax.set(xlim=(0,2),ylim=(0,4),xlabel=r'$\theta/\pi$',ylabel=r'$c_{\rm ent}$',title='Spin 1: bilinear-biquadratic');ax.legend(ncol=3,fontsize=8)
ax=axes[1,2]
for fam,expected in [('uls3',2),('uls4',3),('tb1',1.5),('tb32',1.8)]:
 rs=sorted([r for r in rows if r['family']==fam and np.isfinite(cf(r)) and (not fam.startswith('uls') or r['L']%int(fam[-1])==0)],key=lambda r:r['L'])
 ax.plot([1/r['L'] for r in rs],[cf(r)-expected for r in rs],'o-',label=fam)
ax.axhline(0,color='black',ls='--',lw=1);ax.set(xlabel='1/L',ylabel=r'$c_{\rm ent}-c_{\rm known}$',title='Calibration bias, cut >= 2');ax.legend(fontsize=8)
fig.suptitle('D3 finite-size survey. Entropy excesses are not CFT identifications. Gray regions have exact polarized ground states.',fontsize=12)
fig.savefig('survey.png');fig.savefig('survey.svg');plt.close(fig)

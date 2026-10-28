import os
os.environ['MPLCONFIGDIR']=os.path.join(os.getcwd(),'mpl_cache')
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pair_odd_screen import energy

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
ks=np.linspace(0,np.pi,301)
for delta,col in [(0,'#2166ac'),(.45,'#4d9221'),(.9,'#b35806')]:
 es=[energy(k,-14/9+delta,0,2,10/9,n=100) for k in ks]
 ax[0,0].plot(ks/np.pi,es,color=col,label=r'$\delta='+str(delta)+'$')
ax[0,0].plot(ks/np.pi,-4*np.cos(ks/2),'--',color='.55',lw=1,label='two-atom threshold at K')
ax[0,0].set(xlabel=r'$K/\pi$',ylabel='Two-body energy',ylim=(-4.08,-2.6),title='Molecular minima disappear under detuning')
ax[0,0].legend(frameon=False,fontsize=8)
rows=json.load(open('ed_scaling_delta_0.9.json'))
rows=[r for r in rows if r['L']>=15]+[json.load(open('ed_L24_results.json'))]
Ls=[r['L'] for r in rows]
ax[0,1].plot(Ls,[r['fit']['c_vn'] for r in rows],'o-',label='von Neumann shape fit',color='#2166ac')
ax[0,1].plot(Ls,[r['fit']['c_renyi2'] for r in rows],'s-',label='second-Renyi shape fit',color='#b35806')
ax[0,1].axhline(2,color='.4',ls='--',label='conjectured bound')
ax[0,1].set(xlabel='L',ylabel='Finite-size estimate',title=r'Large shape fits are not a settled $c$')
ax[0,1].legend(frameon=False,fontsize=8)
spec=json.load(open('momentum_scaling_results.json'))
gaps=[]
for r in spec:
 if r['L']<15:continue
 es={q['k']:q['energies'][0] for q in r['levels']}
 gaps.append((r['L'],es[1]-es[0]))
r=json.load(open('ed_L24_results.json'));es={q['k']:q['energies'][0] for q in r['levels']};gaps.append((24,es[1]-es[0]))
x=np.array([z[0] for z in gaps]);y=np.array([z[1] for z in gaps])
ax[1,0].loglog(x,y,'o-',color='#2166ac',label='lowest small-momentum gap')
xx=np.linspace(15,25,100);ref=y[-2]
ax[1,0].loglog(xx,ref*(21/xx),ls='--',color='.4',label=r'$L^{-1}$ guide')
ax[1,0].loglog(xx,ref*(21/xx)**2,ls=':',color='#b35806',label=r'$L^{-2}$ guide')
ax[1,0].set(xlabel='L',ylabel=r'$\Delta E(2\pi/L)$',title='Size drift prevents a clean exponent')
ax[1,0].set_xticks([15,18,21,24],labels=['15','18','21','24'])
from matplotlib.ticker import NullFormatter
ax[1,0].xaxis.set_minor_formatter(NullFormatter())
ax[1,0].legend(frameon=False,fontsize=8)
points=[]
for f in ['periodic_L24_N8_d09.json','periodic_L24_refined.json','periodic_L24_highchi.json']:
 rr=json.load(open(f))
 for r in rr:
  points.append((r['chi'],r['entropy'][11]['S']))
# Last completed sweep at each bond dimension.
vals={chi:S for chi,S in points};chi=sorted(vals)
ax[1,1].plot(chi,[vals[k] for k in chi],'o-',color='#2166ac',label='periodic DMRG, L=24')
ax[1,1].axhline(2.4646300124430667,ls='--',color='.4',label='exact diagonalization')
ax[1,1].set(xlabel='Bond dimension',ylabel='Half-chain von Neumann entropy',title='Finite bond dimension suppresses entropy')
ax[1,1].legend(frameon=False,fontsize=8)
fig.suptitle('C12 qutrit candidate: established two-body structure, inconclusive many-body scaling',fontsize=12)
fig.savefig('diagnostics.png',dpi=180)

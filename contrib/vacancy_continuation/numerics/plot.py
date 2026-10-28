"""Optional figure. Requires matplotlib in addition to numpy."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
result=json.loads((HERE/'summary.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold'})
fig,ax=plt.subplots(1,2,figsize=(10.8,4.4),layout='constrained')
fig.suptitle('Vacancy + XX payload: independent microscopic checks',fontsize=15,fontweight='bold')
colors=['#2d6a9f','#d16b2e','#54814f'];x=np.log(np.array(result['sizes'])/8)
for color,run in zip(colors,result['runs']):
    y=np.array(run['half_chain_entropy']);y-=y[0]
    ax[0].plot(x,y,'o-',color=color,label=f"J = {run['J']:.2f}; c fit = {run['c_entropy_three_sizes']:.3f}",lw=1.4,markersize=5)
ax[0].plot([0,x[-1]],[0,2*x[-1]/3],'k--',label='c = 2 reference',lw=1.3)
ax[0].set(title='Physical half-chain entanglement',xlabel='log(L / 8)',ylabel='S(L/2) − S(4) at L=8')
ax[0].set_xticks(x,[f'log({L}/8)' for L in result['sizes']]);ax[0].legend(frameon=False,fontsize=8,loc='upper left')
js=np.array([r['J'] for r in result['runs']]);response=[r['response_velocity_sum_L16'] for r in result['runs']];casimir=[r['casimir_velocity_sum_Lminus1_Lminus3_fit'] for r in result['runs']]
ax[1].plot(js,response,'o-',color='#2d6a9f',markersize=7,label='Flux + number responses, L=16')
ax[1].plot(js,casimir,'s',mfc='white',mec='#d16b2e',mew=1.7,markersize=6,label='Ground-energy Casimir fit')
grid=np.linspace(0,.22,100);ax[1].plot(grid,2+grid*(2-12/np.pi**2),'k--',lw=1.2,label='First order in J')
ax[1].set(title='Two speeds independently calibrated',xlabel='Payload coupling J / t',ylabel='v_charge + v_spin (t = 1)',xlim=(.025,.22))
ax[1].legend(frameon=False,fontsize=8,loc='upper left')
for a in ax:a.grid(alpha=.15)
fig.supxlabel('Full d=3 model, L=8,12,16. Finite-size evidence for c=2; this saturates c ≤ d−1.',fontsize=9)
fig.savefig(HERE/'vacancy_control.png',dpi=180)
print(HERE/'vacancy_control.png')

"""Summarize only completed finite-chain runs; never infer a CFT from a fit alone."""
import json,math,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/cd-matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent
def get(name,L):
    if name=='c181153':
        paths=[BASE/'results'/f'{name}_L{L}_{parity}.json' for parity in ('even','odd')]
        if not all(p.exists() for p in paths):return None
        return min([json.loads(p.read_text()) for p in paths],key=lambda a:a['E_shifted'])
    p=BASE/'results'/f'{name}_L{L}.json'
    return json.loads(p.read_text()) if p.exists() else None

def exact_xx_entropy(L):
    j=np.arange(1,L+1)[:,None];k=np.arange(1,L//2+1)[None,:]
    # Filling the opposite half of this bipartite band has identical entropy.
    U=np.sqrt(2/(L+1))*np.sin(np.pi*j*k/(L+1));C=U@U.T
    out=[]
    for cut in range(1,L):
        v=np.linalg.eigvalsh(C[:cut,:cut]);v=np.clip(v,1e-15,1-1e-15)
        out.append(float(-np.sum(v*np.log(v)+(1-v)*np.log(1-v))))
    return np.array(out)

records=[];pairs=[];exact=[]
for name in ['c183','dyadic','original','c181153','xx','two_xx']:
    for L in [12,32,64,128]:
        a=get(name,L)
        if a is None:continue
        st=a['sweep_stats'];n=max(1,len(st['E'])-4)
        records.append(dict(name=name,L=L,parity=a.get('parity'),S=a['S_middle3'],density=a['bulk_density'],
                            E=a['E_shifted'],chi=a['chi'],profile_fits=a['profile_fits'],
                            status=('unconverged; nominal parity selection only' if name=='c181153' and L==128 else 'finite-size diagnostic; see convergence history'),
                            final_four_energy_span=max(st['E'][n:])-min(st['E'][n:]),
                            last_mean_entropy_change=st['Delta_S'][-1],last_discarded=st['max_trunc_err'][-1]))
        if name in ('xx','two_xx'):
            if name=='xx':S=exact_xx_entropy(L)
            else:
                s=np.r_[0.,exact_xx_entropy(L//2),0.]
                S=np.array([s[(cut+1)//2]+s[cut//2] for cut in range(1,L)])
            exact.append(dict(name=name,L=L,exact_middle3=float(S[L//2-2:L//2+1].mean()),
                              max_entropy_error=float(np.max(np.abs(S-np.array(a['entropy']))))))
    for L in [32,64]:
        a,b=get(name,L),get(name,2*L)
        if a is not None and b is not None:
            pairs.append(dict(name=name,L1=L,L2=2*L,c_two=6*(b['S_middle3']-a['S_middle3'])/math.log(2),
                              density1=a['bulk_density'],density2=b['bulk_density']))
sector_pairs=[]
for L in [32,64]:
    for p in ('even','odd'):
        for q in ('even','odd'):
            f=BASE/'results'/f'c181153_L{L}_{p}.json';g=BASE/'results'/f'c181153_L{2*L}_{q}.json'
            if f.exists() and g.exists():
                a,b=json.loads(f.read_text()),json.loads(g.read_text())
                sector_pairs.append(dict(L1=L,L2=2*L,parity1=p,parity2=q,c_two=6*(b['S_middle3']-a['S_middle3'])/math.log(2)))
summary=dict(status='No established c>1: parity-family c~2 at 32-64 is not sustained at unconverged L128',records=records,pairs=pairs,exact_controls=exact,parity_sensitivity=sector_pairs)
(BASE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,(ax,bx)=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
colors={'c183':'#2367a3','dyadic':'#d27625','original':'#888888','c181153':'#7b4eaa'}
labels={'c183':'Transverse DM reference','dyadic':'Rational transverse DM','original':'Original V=40 point','c181153':'Parity family, phase 1.53'}
for name in ['c183','dyadic','original','c181153']:
    data=[r for r in records if r['name']==name and r['L']>=32]
    if len(data)<2:continue
    x=np.array([r['L'] for r in data]);y=np.array([r['S'] for r in data]);y-=y[0]
    ax.plot(x,y,'-',color=colors[name],label=labels[name])
    ax.scatter(x[x<128],y[x<128],color=colors[name],s=30)
    ax.scatter(x[x==128],y[x==128],facecolors='white',edgecolors=colors[name],s=40,zorder=3)
grid=np.linspace(32,128,200)
for c in [1,2]:ax.plot(grid,c*np.log(grid/32)/6,'--',color='#666666',alpha=.7,label=f'c={c} slope')
ax.set_xscale('log',base=2);ax.set_xticks([32,64,128],['32','64','128'])
ax.set_xlabel('Sites L');ax.set_ylabel('Change in central entropy from L=32')
ax.set_title('Full microscopic chains: size scaling');ax.legend(fontsize=8,loc='lower right')
ax.text(.02,.98,'Hollow: L=128, sweep limit reached',transform=ax.transAxes,ha='left',va='top',fontsize=8)
short=[p for p in pairs if p['name'] in ['xx','two_xx','c183','dyadic','c181153']]
names=[f"{ {'xx':'XX control','two_xx':'Two XX legs','c183':'Transverse DM','dyadic':'Rational DM','c181153':'Parity family'}[p['name']]}  {p['L1']}→{p['L2']}" for p in short]
vals=[p['c_two'] for p in short]
bx.barh(names,vals,color=[colors.get(p['name'],'#9da5af') for p in short])
bx.axvline(1,color='#444444',ls='--',lw=1);bx.axvline(2,color='#444444',ls=':',lw=1)
for i,v in enumerate(vals):bx.text(v+.025,i,f'{v:.2f}',va='center',fontsize=10)
bx.set_xlim(0,max(2.4,max(vals,default=2)+.3));bx.set_xlabel('Two-size entropy slope (not extrapolated c)')
bx.set_title('The parity-family signal falls at L=128')
fig.savefig(BASE/'scaling.png',dpi=180)
fig.savefig(BASE/'scaling.pdf')
print(json.dumps(pairs,indent=2))

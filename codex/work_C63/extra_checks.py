import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from ed_staggered import low,sector,entropy_profile

# Check every conserved-charge sector at the largest main-scan size.
points=[(0.,0.),(.5,0.),(.75,.5),(.9,.5),(1.3,2.),(1.35,2.),(1.4,2.),(2.2,4.),(2.22,4.)]
L=18;out=[]
for u,b in points:
    energies=[]
    for n in range(L//2+1):
        H,_=sector(L,n,u,b);e,_,r=low(H,1);energies.append(float(e[0]))
    item=dict(u=u,b=b,L=L,energies_N_0_to_half=energies,minimizer=int(np.argmin(energies)))
    assert item['minimizer']==L//2,item
    out.append(item);print(json.dumps(item),flush=True)
Path('all_sector_checks.json').write_text(json.dumps(out,indent=2)+'\n')
# Independent c=1/2 and additive c=3/2 entropy controls.
controls=[]
for L in [8,10,12,14,16,18]:
    st=np.arange(1<<L,dtype=np.int64);st=st[np.bitwise_count(st)%2==0]
    z=1-2*((st[:,None]>>np.arange(L))&1)
    rr=[np.arange(len(st))];cc=[rr[0]];dd=[-np.sum(z,axis=1).astype(float)]
    for j in range(L):
        rr.append(np.searchsorted(st,st^((1<<j)|(1<<((j+1)%L)))));cc.append(cc[0]);dd.append(-np.ones(len(st)))
    H=coo_matrix((np.concatenate(dd),(np.concatenate(rr),np.concatenate(cc))),shape=(len(st),len(st))).tocsr()
    e,v,r=low(H,1);S,cS,err=entropy_profile(v[:,0],st,L)
    Hx,sx=sector(L,L//2,0,0);ex,vx,rx=low(Hx,1);Sx,cx,_=entropy_profile(vx[:,0],sx,L)
    item=dict(L=L,E_Ising=float(e[0]),exact_energy_error=float(abs(e[0]+2/np.sin(np.pi/(2*L)))),S_Ising=S,c_Ising=cS,S_product=(np.array(S)+Sx).tolist(),c_product=cS+cx,E_product=float(e[0]+ex[0]),residual=max(r,rx))
    controls.append(item);print(json.dumps({k:item[k] for k in ['L','c_Ising','c_product','exact_energy_error','residual']}),flush=True)
Path('controls.json').write_text(json.dumps(controls,indent=2)+'\n')

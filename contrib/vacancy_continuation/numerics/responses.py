"""Independent finite-size stiffness and compressibility of both U(1)s."""
import argparse,json,time
import numpy as np
from scipy.sparse import coo_matrix
from xx_vacancy import HERE,build,ground

def flux_hamiltonian(codes,L,J,phi,charges):
    powers=3**np.arange(L,dtype=np.int64)
    digits=(codes[:,None]//powers[None,:])%3
    rows=[];cols=[];data=[];q=np.array(charges)
    for x in range(L):
        y=(x+1)%L;a=digits[:,x];b=digits[:,y]
        for coupling,mask in ((1.,(a==0)^(b==0)),(J,(a+b==3)&(a!=0)&(b!=0))):
            ii=np.flatnonzero(mask)
            target=codes[ii]+(b[ii]-a[ii])*(powers[x]-powers[y])
            jj=np.searchsorted(codes,target)
            rows.append(jj);cols.append(ii)
            data.append(-coupling*np.exp(1j*phi*(q[a[ii]]-q[b[ii]])/L))
    H=coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(len(codes),len(codes))).tocsr()
    delta=H-H.conj().T
    assert delta.nnz==0 or np.max(abs(delta.data))<1e-13
    return H

def run(L,js):
    started=time.time(); n=L//4
    codes,mats=build(L,n,n)
    out={'L':L,'dimension':len(codes),'runs':[]}
    for J in js:
        H=-mats[0]-J*mats[1]
        ev,vs,res=ground(H);e0=ev[0];psi=vs[:,0]
        item={'J':J,'E0':float(e0),'ground_residual':res[0],'flux':{}}
        for name,charges in [('charge',[0,1,1]),('spin',[0,-.5,.5])]:
            vals=[]
            for phi in [.04,.02]:
                HF=flux_hamiltonian(codes,L,J,phi,charges)
                e,v,r=ground(HF,psi.astype(complex))
                vals.append({'phi':phi,'E':float(e[0]),'residual':r[0],
                             'vK':float(2*np.pi*L*(e[0]-e0)/phi**2)})
                del HF
            item['flux'][name]=vals
        energies={}
        for name,u,d in [('spin_plus',n+1,n-1),('charge_plus',n+1,n+1),('charge_minus',n-1,n-1)]:
            cs,ms=build(L,u,d)
            e,v,r=ground(-ms[0]-J*ms[1]);energies[name]={'energy':float(e[0]),'dimension':len(cs),'residual':r[0]}
            del cs,ms,v
        item['sectors']=energies
        vk_spin=2*L*(energies['spin_plus']['energy']-e0)/np.pi
        vk_charge=L*(energies['charge_plus']['energy']+energies['charge_minus']['energy']-2*e0)/(4*np.pi)
        for name,v_over_k in [('spin',vk_spin),('charge',vk_charge)]:
            vals=item['flux'][name]
            v_times_k=(4*vals[1]['vK']-vals[0]['vK'])/3
            assert v_times_k>0 and v_over_k>0
            item[name]={'v_over_K':float(v_over_k),'v_times_K_phi_extrap':float(v_times_k),
                        'v':float(np.sqrt(v_times_k*v_over_k)),'K':float(np.sqrt(v_times_k/v_over_k))}
        item['velocity_sum']=item['spin']['v']+item['charge']['v']
        item['elapsed_seconds']=time.time()-started
        out['runs'].append(item)
        (HERE/f'response_L{L}.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({'L':L,'J':J,'charge':item['charge'],'spin':item['spin'],'seconds':time.time()-started}),flush=True)
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--L',type=int,default=12);ap.add_argument('--J',type=float,nargs='+',default=[.2,.1,.05]);a=ap.parse_args()
    run(a.L,a.J)

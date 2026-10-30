"""Independent finite-dimensional temporal-transfer identities."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh
from checks import dense_xx,dense_histories

ROOT=Path(__file__).resolve().parent

def diagonalized_transfer(ev,vec,tau,mask):
    e=(vec*np.exp(-tau*ev/2))@vec.T
    keep=np.flatnonzero(mask.ravel())
    full=np.kron(e,e)
    t=full[np.ix_(keep,keep)]
    u=e.ravel()[keep]
    vals,basis=eigh(t)
    assert vals[0]>-1e-10
    weights=(basis.T@u)**2
    top=vals>=vals[-1]*(1-1e-11)
    return vals,weights,float(weights[top].sum()),float(u@u)

def main():
    ev,v=dense_xx(2)
    evd,vd=dense_xx(4)
    h=(v*ev)@v.T
    h0=np.kron(h,np.eye(4))+np.kron(np.eye(4),h)
    e0,v0=eigh(h0)
    # Tensor labels are (a,b), each ranging over the four states of a 2-ring.
    pair=np.arange(16)
    ca=(pair//4)%2; cb=pair%2
    maskS=(ca[:,None]==ca[None,:])&(ca[:,None]==cb[:,None])&(ca[None,:]==cb[None,:])
    labelD=(pair%2)+2*((pair//4)%2)
    maskD=labelD[:,None]==labelD[None,:]
    rows=[]
    for tau in [.2,1.,3.]:
        vs,ws,topS,normS=diagonalized_transfer(e0,v0,tau,maskS)
        vd0,wd,topD,normD=diagonalized_transfer(evd,vd,tau,maskD)
        rate=float(np.log(vs[-1]/vd0[-1]))
        boundary=float(np.log(normS/topD))
        for N in [1,2,3,4]:
            beta=N*tau
            q=dense_histories(ev,v,beta,N,[0])*np.exp(-beta*ev[0])
            qd=dense_histories(evd,vd,beta,N,[0,2])*np.exp(-beta*evd[0])
            s=float(ws@(vs**(N-1))); d=float(wd@(vd0**(N-1)))
            err=max(abs(s-q@q)/s,abs(d-qd.sum())/d)
            assert err<1e-11,err
            g=np.log(s/d); upper=boundary+(N-1)*rate
            assert g<=upper+1e-10
            rows.append(dict(tau=tau,N=N,error=err,G=float(g),upper=float(upper),
                             rate=rate,boundary=boundary))
    (ROOT/'transfer_results.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('transfer checks',len(rows),'max error',max(r['error'] for r in rows))

if __name__=='__main__': main()

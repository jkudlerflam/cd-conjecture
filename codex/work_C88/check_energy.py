"""Resolved color self-energy, including nonzero payload bits."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,itertools
import numpy as np
from scipy.sparse.linalg import eigsh
from check_lift import closure,run
from bond13 import params,incoming,outgoing

def main():
    out=[]
    for bits in itertools.product((0,1),repeat=2):
        seed=[('R',0)]+[('R',a) for a in bits]+[('e',0),('V',0)]
        # Color the vacancy after E by state 0.
        H,pen,bad,states=closure(seed,params(n=3,K=0.),0.,True)
        ids=np.flatnonzero(pen==0);P=H[ids][:,ids]
        es,vs=eigsh(P,k=1,which='SA',tol=1e-13,v0=np.ones(len(ids)))
        e=es[0];v=vs[:,0]
        ct=[];cv=[]
        for j in ids:
            w=states[j];L=len(w)
            ct.append(sum(len(incoming(a))-1 for i,a in enumerate(w) if a[0]!='V' and w[(i+1)%L][0]=='V'))
            cv.append(sum(outgoing(w[(i-1)%L])==1 and a==('R',0) and w[(i+1)%L][0] in ('R','e') for i,a in enumerate(w)))
        et=float(np.vdot(v,np.array(ct)*v).real)
        ew=float(np.vdot(v,np.array(cv)*v).real)
        records=[]
        for K in (100.,1000.,10000.):
            rec=run(H,pen,bad,K)
            rec['K_energy_shift']=float(K*(rec['energy']-e));records.append(rec)
        assert abs(records[-1]['K_energy_shift']+et+ew)<.003
        out.append(dict(bits=list(bits),dimension=len(pen),legal_dimension=len(ids),
          constrained_energy=float(e),color_shift_coefficient=-et,
          vacuum_shift_coefficient=-ew,runs=records))
        print(json.dumps(out[-1]),flush=True)
    H,pen,bad,states=closure([('L',1),('R',0),('e',0),('V',0)],params(n=2,K=0.),.1,True,J=.01,t=1.)
    lowJ=[]
    for K in (10000.,100000.,1000000.):
        r=run(H,pen,bad,.01*K);r['K']=K;r['kappa']=.01*K;r['J']=.01;lowJ.append(r)
    json.dump(dict(energy=out,small_J=lowJ),open('energy_checks.json','w'),indent=2)
    print(json.dumps(dict(small_J=lowJ)),flush=True)
if __name__=='__main__':main()

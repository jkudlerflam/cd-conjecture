import numpy as np,signal
from periodic_dmrg import PeriodicDMRG
from simple_dmrg import qutrit_model

def random_charge_seed(obj,N,seed=912):
    rng=np.random.default_rng(seed);L=obj.L;d=obj.d
    obj.q=[np.arange(max(0,N-(d-1)*(L-i)),min((d-1)*i,N)+1) for i in range(L+1)]
    obj.mps=[]
    for i in range(L):
        ql,qr=obj.q[i],obj.q[i+1]
        M=np.zeros((len(ql),d,len(qr)))
        for a,x in enumerate(ql):
            for s in range(d):
                where=np.flatnonzero(qr==x+s)
                if len(where):M[a,s,where[0]]=rng.uniform(.8,1.2)
        obj.mps.append(M)
    # One vector per charge makes distinct rows orthogonal automatically.
    for i in range(L-1,0,-1):
        norms=np.sqrt(np.sum(obj.mps[i]**2,axis=(1,2)))
        assert np.all(norms>0)
        obj.mps[i]/=norms[:,None,None]
        obj.mps[i-1]*=norms[None,None,:]
    obj.mps[0]/=np.linalg.norm(obj.mps[0])
    obj.rebuild_right()

if __name__=='__main__':
    signal.alarm(2100)
    h,o=qutrit_model(.9,mu=0)
    obj=PeriodicDMRG(48,16,h,o)
    random_charge_seed(obj,16)
    obj.run([24,32,48,64,96,128,128,128,128,128],output='periodic_L48_N16_d09_random',tol=1e-9)

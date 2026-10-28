"""Translation-orbit reduction, with explicit comparison against full ED."""
import numpy as np,json,time
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from qutrit_ed import basis,local_h

def rotate(states,L):
    if states.dtype==object:
        return ((states<<2)&((1<<(2*L))-1))|(states>>(2*(L-1)))
    mask=np.uint64((1<<(2*L))-1)
    return ((states<<np.uint64(2))&mask)|(states>>np.uint64(2*(L-1)))

def canonical(states,L):
    best=states.copy();shift=np.zeros(len(states),np.int16);cur=states.copy()
    for j in range(1,L):
        cur=rotate(cur,L)
        smaller=cur<best;best[smaller]=cur[smaller];shift[smaller]=j
    return best,shift

class OrbitHamiltonian:
    def __init__(self,L,N,custom_basis=None,custom_h=None,**pars):
        self.L=L;self.N=N;start=time.time()
        self.custom_basis=custom_basis
        states=basis(L,N) if custom_basis is None else custom_basis
        canon,_=canonical(states,L)
        self.reps=states[canon==states]
        self.full_dim=len(states)
        del states,canon
        reps=self.reps;D=len(reps)
        cur=reps.copy();period=np.full(D,L,dtype=np.int16)
        for j in range(1,L):
            cur=rotate(cur,L)
            hit=(cur==reps)&(period==L);period[hit]=j
        self.period=period
        h=local_h(**pars) if custom_h is None else custom_h;diag=np.zeros(D);rows=[];cols=[];vals=[];powers=[]
        for i in range(L):
            j=(i+1)%L
            a=((reps>>(2*i))&3).astype(np.int8)
            b=((reps>>(2*j))&3).astype(np.int8)
            loc=3*a+b;diag+=np.diag(h)[loc]
            for old in range(9):
                for new in range(9):
                    value=h[new,old]
                    if new==old or value==0:continue
                    ix=np.flatnonzero(loc==old)
                    if not len(ix):continue
                    aa,bb=divmod(new,3);ao,bo=divmod(old,3)
                    delta=(aa-ao)*(1<<(2*i))+(bb-bo)*(1<<(2*j))
                    if reps.dtype==object:target=reps[ix]+delta
                    else:target=reps[ix]+np.uint64(delta) if delta>=0 else reps[ix]-np.uint64(-delta)
                    can,shift=canonical(target,L);dest=np.searchsorted(reps,can)
                    assert np.all(reps[dest]==can)
                    rows.append(dest);cols.append(ix);powers.append(shift)
                    vals.append(value*np.sqrt(period[ix]/period[dest]))
        rows.append(np.arange(D));cols.append(np.arange(D));powers.append(np.zeros(D,dtype=np.int16));vals.append(diag)
        self.rows=np.concatenate(rows);self.cols=np.concatenate(cols);self.powers=np.concatenate(powers);self.vals=np.concatenate(vals)
        self.build_seconds=time.time()-start

    def sector(self,k):
        valid=(self.period*k)%self.L==0
        remap=np.cumsum(valid)-1
        keep=valid[self.rows]&valid[self.cols]
        vals=self.vals[keep]*np.exp(-2j*np.pi*k*self.powers[keep]/self.L)
        D=sum(valid)
        H=coo_matrix((vals,(remap[self.rows[keep]],remap[self.cols[keep]])),shape=(D,D)).tocsr()
        return H,valid

    def solve(self,k,nev=3,tol=1e-10):
        H,valid=self.sector(k);start=time.time()
        err=H-H.T.conj()
        assert not err.nnz or np.max(abs(err.data))<1e-10
        if H.shape[0]<max(nev+2,12):
            es,vs=np.linalg.eigh(H.toarray());es,vs=es[:nev],vs[:,:nev]
        else:
            rng=np.random.default_rng(100+k)
            es,vs=eigsh(H,k=nev,which='SA',tol=tol,v0=rng.normal(size=H.shape[0])+1j*rng.normal(size=H.shape[0]),ncv=max(2*nev+1,24))
            ix=np.argsort(es);es,vs=es[ix],vs[:,ix]
        return es,vs,dict(k=k,dim=H.shape[0],energies=es.tolist(),seconds=time.time()-start,residual=float(np.linalg.norm(H@vs[:,0]-es[0]*vs[:,0])))

    def reconstruct(self,k,v):
        states=basis(self.L,self.N) if self.custom_basis is None else self.custom_basis;can,shift=canonical(states,self.L)
        orb=np.searchsorted(self.reps,can)
        valid=(self.period*k)%self.L==0;remap=np.cumsum(valid)-1
        good=valid[orb]
        full=np.zeros(len(states),complex)
        full[good]=v[remap[orb[good]]]*np.exp(2j*np.pi*k*shift[good]/self.L)/np.sqrt(self.period[orb[good]].astype(float))
        return states,full

if __name__=='__main__':
    from qutrit_ed import matrix
    L,N=8,4
    ob=OrbitHamiltonian(L,N,U=-14/9+.9)
    vals=[]
    for k in range(L):
        es,vs,r=ob.solve(k,nev=3);vals+=es.tolist()
        states,v=ob.reconstruct(k,vs[:,0]);H,_=matrix(L,N,U=-14/9+.9)
        residual=np.linalg.norm(H@v-es[0]*v)
        print(k,es,residual)
        assert residual<1e-8
    E=eigsh(H,k=1,which='SA')[0][0]
    print('Ground full / all momentum sectors',E,min(vals))
    assert abs(E-min(vals))<1e-8

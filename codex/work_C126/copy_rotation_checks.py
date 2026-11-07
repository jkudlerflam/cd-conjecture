"""Independent small matrix checks for the C126 derivations."""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, expm_multiply
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent


def local_sum(h,n,edges=None):
    edges=edges or [(x,(x+1)%n) for x in range(n)]
    states=np.arange(2**n,dtype=np.int64)
    rows=[];cols=[];values=[]
    for x,y in edges:
        px=1<<(n-1-x);py=1<<(n-1-y)
        a=(states//px)%2;b=(states//py)%2
        old=2*a+b
        for out in range(4):
            value=h[out,old]
            mask=abs(value)>1e-14
            target=states[mask]+(out//2-a[mask])*px+(out%2-b[mask])*py
            rows.append(target);cols.append(states[mask]);values.append(value[mask])
    return sparse.coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(2**n,2**n)).tocsr()


def ground(h,n):
    H=local_sum(h,n)
    e,v=eigsh(H,k=1,which='SA',tol=1e-12,v0=np.cos(np.arange(2**n)*.719+.4))
    return float(e[0]),v[:,0],float(np.linalg.norm(H@v[:,0]-e[0]*v[:,0]))


def identity_check():
    rng=np.random.default_rng(126)
    X=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));h=(X+X.conj().T)/2
    errors=[]
    for L,N in [(2,2),(3,2),(2,3)]:
        # Rung-major coordinates on the twisted system.
        edges=[(x*N+r,(x+1)*N+r) for x in range(L-1) for r in range(N)]
        edges += [((L-1)*N+r,(r+1)%N) for r in range(N)]
        a=np.linalg.eigvalsh(local_sum(h,L*N,edges).toarray())
        b=np.linalg.eigvalsh(local_sum(h,L*N).toarray())
        errors.append(float(max(abs(a-b))))
    return errors


def free_checks():
    sx=np.array([[0,1],[1,0]])/2
    sy=np.array([[0,-1j],[1j,0]])/2
    sz=np.diag([.5,-.5]);I=np.eye(2)
    rows=[]
    for jx,jy,g in [(1,1,0),(1,1,.7),(1,.3,.7),(1,0,.5),(1,-.4,.23)]:
        hop=(jx*np.kron(sx,sx)+jy*np.kron(sy,sy)).real
        h=hop+g/2*(np.kron(sz,I)+np.kron(I,sz))
        ev=np.linalg.eigvalsh(h);M=float((ev[-1]-ev[0])/2)
        for L,N in [(2,2),(4,2),(6,2),(4,3)]:
            e,v,res=ground(h,L);en,_,resn=ground(h,L*N)
            ehop=float(v@(local_sum(hop,L)@v)/L)
            phases=2*np.pi*(np.arange(N)-(N-1)/2)/N
            trial=L*ehop*(np.sum(np.cos(phases/L))-N)
            cap=np.pi**2*M/(6*L)*(N-1/N)
            assert en-N*e <= trial+1e-9 and trial<=cap+1e-9
            if jx==jy and g==0: assert abs(trial-(en-N*e))<1e-9
            rows.append(dict(jx=jx,jy=jy,g=g,L=L,N=N,exact=en-N*e,trial=trial,
                             norm_bound=cap,M=M,residual=max(res,resn)))
    return rows


def cat_direct():
    # Direct full-vector evaluation, independent of the doubly stochastic formula.
    rng=np.random.default_rng(222)
    N=2;L=3;D=2**N
    state=np.zeros([D]*L,dtype=complex)
    for a in range(D): state[(a,)*L]=1/np.sqrt(D)
    units=[]
    for x in range(L):
        z=rng.normal(size=(D,D))+1j*rng.normal(size=(D,D))
        u,_=np.linalg.qr(z);units.append(u)
        state=np.tensordot(u,state,axes=(1,x));state=np.moveaxis(state,0,x)
    Z=np.diag([1.,-1.]);h=-np.kron(Z,Z)
    edges=[(x*N+r,(x+1)*N+r) for x in range(L-1) for r in range(N)]
    edges += [((L-1)*N+r,(r+1)%N) for r in range(N)]
    v=state.ravel();direct=float((v.conj()@(local_sum(h,L*N,edges)@v)).real+N*L)
    zz=1-2*((np.arange(D)[:,None]//(2**np.arange(N-1,-1,-1)))%2)
    bars=[abs(u).T**2@zz for u in units]
    formula=sum(N-np.sum(bars[x]*bars[x+1])/D for x in range(L-1))
    formula+=N-np.sum(bars[-1]*np.roll(bars[0],1,axis=1))/D
    assert abs(direct-formula)<1e-11
    return dict(direct=direct,formula=float(formula),error=abs(direct-formula))


def explicit_fermion_trial():
    sx=np.array([[0,1],[1,0]])/2
    sy=np.array([[0,-1j],[1j,0]])/2
    sz=np.diag([.5,-.5]);I=np.eye(2);L=4;n=2*L
    rows=[]
    for jy,g in [(1.,0.),(1.,.7),(.3,.7),(0.,.5),(-.4,.23)]:
        hop=(np.kron(sx,sx)+jy*np.kron(sy,sy)).real
        h=hop+g/2*(np.kron(sz,I)+np.kron(I,sz))
        e,psi,_=ground(h,L)
        parity=np.array([(-1)**int(i).bit_count() for i in range(2**L)])
        p=1 if np.dot(psi**2,parity)>0 else -1
        psi=psi*(parity==p);psi/=np.linalg.norm(psi)
        assert np.linalg.norm(local_sum(h,L)@psi-e*psi)<1e-10
        eta=-p
        G=np.zeros((2**n,2**n))
        for state in range(2**n):
            digits=[(state>>(n-1-i))&1 for i in range(n)]
            for x in range(L):
                a=x;b=L+x;angle=x*np.pi/(2*L)
                if digits[a]==digits[b]:continue
                target=state^(1<<(n-1-a))^(1<<(n-1-b))
                sign=(-1)**sum(digits[a+1:b])
                G[target,state]+=angle*sign*(1 if digits[b] else -1)
        omega=np.kron(psi,psi)
        trial=expm_multiply(sparse.csr_matrix(G),omega)
        gauge=np.array([eta**int(i%(2**L)).bit_count() for i in range(2**n)])
        trial*=gauge
        excess=float(trial@(local_sum(h,n)@trial)-2*e)
        ehop=float(psi@(local_sum(hop,L)@psi)/L)
        formula=2*L*ehop*(np.cos(np.pi/(2*L))-1)
        assert abs(excess-formula)<1e-10,(jy,g,excess,formula)
        rows.append(dict(jy=jy,g=g,L=L,parity=p,direct=excess,formula=formula,error=abs(excess-formula)))
    return rows


def cat_all_N():
    records=[]
    rng=np.random.default_rng(126)
    for N in range(2,7):
        D=2**N;P=np.zeros((D,D))
        codes=np.arange(D)
        targets=((codes<<1)&(D-1))+(codes>>(N-1))
        P[targets,codes]=1
        z=1-2*((codes[:,None]//(2**np.arange(N-1,-1,-1)))%2)
        phases=2*np.pi*np.arange(N)/N
        phases=np.where(phases>np.pi,phases-2*np.pi,phases)
        powers=[np.linalg.matrix_power(P,k) for k in range(N)]
        def fractional(t):
            return sum(np.mean(np.exp(1j*phases*(t-k)))*powers[k] for k in range(N))
        def cost(units):
            bars=[abs(u).T**2@z for u in units]
            total=sum(N-np.sum(bars[x]*bars[x+1])/D for x in range(len(units)-1))
            # Direction is immaterial here, provided it matches P at the seam.
            total+=N-np.sum(bars[-1]*np.roll(bars[0],-1,axis=1))/D
            return float(total)
        random_min=float('inf')
        for _ in range(3):
            units=[np.linalg.qr(rng.normal(size=(D,D))+1j*rng.normal(size=(D,D)))[0] for x in range(4)]
            random_min=min(random_min,cost(units))
        identity=cost([np.eye(D)]*4)
        linear=cost([fractional(x/8) for x in range(8)])
        uniform=cost([fractional(.5)]*8)
        assert abs(identity-N)<1e-10 and random_min>=N-1e-10
        records.append(dict(N=N,identity=identity,random_min=random_min,linear_L8=linear,
                            uniform_L8=uniform,linear_per_copy=linear/N))
    return records


def stabilizer_check():
    pauli=[np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.])]
    h=sum(np.kron(p,p) for p in pauli[1:]).real/4
    L=6;e,v,_=ground(h,L);m=v.reshape(4,-1);rho=m@m.T
    T=np.array([[np.trace(rho@np.kron(a,b)).real for b in pauli] for a in pauli])
    r=4*e/(3*L)
    error=float(np.max(abs(T-np.diag([1.,r,r,r]))))
    assert error<1e-10
    return dict(L=L,r=r,correlation_error=error,eigenvalues=np.linalg.eigvalsh(np.kron(T,T)).tolist())


def multimode_trial_check():
    rng=np.random.default_rng(321)
    I=np.eye(2);X=np.array([[0,1],[1,0]]);Y=np.array([[0,-1j],[1j,0]]);Z=np.diag([1.,-1.])
    gamma=[np.kron(X,I),np.kron(Y,I),np.kron(Z,X),np.kron(Z,Y)]
    P=np.kron(Z,Z);modes=2;d=4
    def ring(h,n):
        states=np.arange(d**n,dtype=np.int64);powers=d**np.arange(n-1,-1,-1,dtype=np.int64)
        rows=[];cols=[];vals=[]
        for x in range(n):
            y=(x+1)%n;a=states//powers[x]%d;b=states//powers[y]%d
            for out in range(d*d):
                vv=h[out,d*a+b];mask=abs(vv)>1e-14
                rows.append(states[mask]+(out//d-a[mask])*powers[x]+(out%d-b[mask])*powers[y])
                cols.append(states[mask]);vals.append(vv[mask])
        return sparse.coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(d**n,d**n)).tocsr()
    records=[]
    for L in [2,3]:
        for sample in range(4):
            B=rng.normal(size=(4,4));A=rng.normal(size=(4,4));A=A-A.T
            onsite=sum(1j*A[a,b]/4*(gamma[a]@gamma[b]) for a in range(4) for b in range(4))
            hop=sum(1j*B[a,b]/2*np.kron(gamma[a]@P,gamma[b]) for a in range(4) for b in range(4))
            h=hop+(np.kron(onsite,np.eye(d))+np.kron(np.eye(d),onsite))/2
            ev,vec=np.linalg.eigh(ring(h,L).toarray());e=float(ev[0]);psi=vec[:,0]
            pp=np.array([(-1)**int(i).bit_count() for i in range(d**L)])
            parity=1 if np.dot(abs(psi)**2,pp)>0 else -1
            psi*=pp==parity;psi/=np.linalg.norm(psi);eta=-parity
            nmode=2*L*modes;states=np.arange(2**nmode,dtype=np.int64)
            rows=[];cols=[];vals=[]
            for x in range(L):
                for orbital in range(modes):
                    a=modes*x+orbital;b=modes*(L+x)+orbital
                    ba=(states>>(nmode-1-a))&1;bb=(states>>(nmode-1-b))&1;mask=ba!=bb
                    selected=states[mask]
                    between=sum((selected>>(nmode-1-j))&1 for j in range(a+1,b))
                    rows.append(selected^(1<<(nmode-1-a))^(1<<(nmode-1-b)));cols.append(selected)
                    vals.append(x*np.pi/(2*L)*((-1.)**between)*(2*bb[mask]-1))
            G=sparse.coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(len(states),len(states))).tocsr()
            trial=expm_multiply(G,np.kron(psi,psi))
            trial*=np.array([eta**int(i%(d**L)).bit_count() for i in states])
            direct=float(np.vdot(trial,ring(h,2*L)@trial).real-2*e)
            ehop=float(np.vdot(psi,ring(hop,L)@psi).real/L)
            formula=2*L*ehop*(np.cos(np.pi/(2*L))-1)
            assert abs(direct-formula)<1e-9,(L,sample,direct,formula)
            records.append(dict(L=L,modes=modes,sample=sample,parity=parity,direct=direct,formula=formula,error=abs(direct-formula)))
    return records


if __name__=='__main__':
    result=dict(status='numerical',seam_spectrum_errors=identity_check(),free=free_checks(),cat=cat_direct(),explicit_free_trial=explicit_fermion_trial(),cat_all_N=cat_all_N(),stabilizer=stabilizer_check(),multimode_trial=multimode_trial_check())
    (ROOT/'copy_rotation_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

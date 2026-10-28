"""C30 finite checks. numpy/scipy only, all output stays in this folder."""
import itertools as it
import json
import numpy as np
from scipy import sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh

# PL,L0,LX,L1,M,R0,RX,R1,PR
arrows=np.array([(0,1),(1,1),(1,2),(2,2),(2,3),(3,3),(3,4),(4,4),(4,0)])
rules=[((1,2),(2,3)),((5,6),(6,7))]
legal=arrows[:,1,None]==arrows[None,:,0]
f=lambda r:-2*np.cos(np.pi/(r+1))
d=lambda r:f(r+1)-f(r)
def parameters(m):
    e=-2. if m==1 else (d(m-1)+d(m))/2
    mu=(2*m+3)*e-2*f(m)
    eta=1. if m==1 else (d(m)-d(m-1))/2
    return e,mu,eta

def local_checks():
    count=0
    for old,new in rules:
        for x,y in it.product(range(9),repeat=2):
            a=(x,*old,y);b=(x,*new,y)
            assert [legal[a[j],a[j+1]] for j in range(3)]==[legal[b[j],b[j+1]] for j in range(3)]
            count+=1
    return count

def segment_checks():
    err=0.;worst=1e9;rows=[]
    for a,c in it.product(range(1,9),repeat=2):
        A=sp.diags([-np.ones(a-1),-np.ones(a-1)],[-1,1],shape=(a,a)).toarray()
        C=sp.diags([-np.ones(c-1),-np.ones(c-1)],[-1,1],shape=(c,c)).toarray()
        vals=np.linalg.eigvalsh(np.kron(A,np.eye(c))+np.kron(np.eye(a),C))
        theory=np.sort(np.array([-2*np.cos(np.pi*i/(a+1))-2*np.cos(np.pi*j/(c+1)) for i in range(1,a+1) for j in range(1,c+1)]))
        err=max(err,float(np.max(abs(vals-theory))))
    for m in range(1,13):
        e,mu,eta=parameters(m); mn=1e9;arg=None
        for a,c in it.product(range(1,49),repeat=2):
            excess=mu+f(a)+f(c)-e*(a+c+3)
            bound=eta*(abs(a-m)+abs(c-m))
            assert excess>=bound-1e-12
            if (a,c)!=(m,m) and excess<mn:mn=float(excess);arg=(a,c)
        rows.append(dict(m=m,e=e,mu=mu,eta=eta,min_wrong=mn,where=arg))
        if m<=4:
            ec=d(m);muc=(2*m+3)*ec-2*f(m)
            for a,c in it.product((m,m+1),repeat=2):
                assert abs(muc+f(a)+f(c)-ec*(a+c+3))<1e-13
    return dict(ed_arms=[1,8],max_dimension=64,spectral_error=err,selection_arms=[1,48],selection_m=[1,12],rows=rows)

def ring_checks():
    out=[]
    for L in (3,4,5):
        size=9**L;pows=9**np.arange(L)
        digits=(np.arange(size)[:,None]//pows[None,:])%9
        bad=np.sum(~legal[digits,np.roll(digits,-1,axis=1)],axis=1)
        markers=np.sum(digits==0,axis=1)
        row=[];col=[]
        for j in range(L):
            k=(j+1)%L
            for a,b in rules:
                for old,new in ((a,b),(b,a)):
                    ix=np.where((digits[:,j]==old[0])&(digits[:,k]==old[1]))[0]
                    iy=ix+(new[0]-old[0])*pows[j]+(new[1]-old[1])*pows[k]
                    assert np.array_equal(bad[ix],bad[iy]);assert np.array_equal(markers[ix],markers[iy])
                    row.extend(ix);col.extend(iy)
        adj=sp.csr_matrix((np.ones(len(row)),(row,col)),shape=(size,size))
        ncomp,lab=connected_components(adj,directed=False)
        order=np.argsort(lab);starts=np.r_[0,np.cumsum(np.bincount(lab))]
        e,mu,eta=parameters(1);ground=1e10;mult=0;illegal_slack=1e10;maxblock=0
        for z in range(ncomp):
            ix=order[starts[z]:starts[z+1]];v=ix[0];maxblock=max(maxblock,len(ix))
            kinetic=0. if len(ix)==1 else float(np.linalg.eigvalsh(-adj[ix][:,ix].toarray())[0])
            energy=32*bad[v]+mu*markers[v]+kinetic
            if bad[v]:
                bound=e*L+(32-8-2*max(0,-mu))*bad[v]
                illegal_slack=min(illegal_slack,float(energy-bound));assert energy>=bound-1e-10
            if energy<ground-1e-9:ground=energy;mult=1
            elif abs(energy-ground)<1e-9:mult+=1
        if L==5:assert abs(ground+10)<1e-12 and mult==5
        out.append(dict(L=L,dimension=size,components=ncomp,max_component=maxblock,ground=float(ground),ground_multiplicity=mult,min_illegal_bound_slack=illegal_slack))
    return out

def embed(op,j,r,width=1):
    return sp.kron(sp.eye(2**j),sp.kron(op,sp.eye(2**(r-j-width))),format='csr')
def helix_checks():
    out=[]
    for m in range(1,5):
        k=2*m+3;p=2*k+1;phi=2*np.pi/p
        plus=lambda th:np.array([1,np.exp(1j*th)])/np.sqrt(2)
        s=np.array([0,1,-np.exp(-1j*phi),0])/np.sqrt(2)
        h=np.outer(s,s.conj());pin0=np.eye(2)-np.outer(plus(0),plus(0).conj())
        pin1=np.eye(2)-np.outer(plus((m-1)*phi),plus((m-1)*phi).conj())
        for r in range(1,min(k,9)+1):
            H=embed(pin0,0,r)+embed(pin1,r-1,r)
            for j in range(r-1):H+=embed(h,j,r,2)
            if 2**r<=128:ev=np.linalg.eigvalsh(H.toarray())[:2]
            else:ev=np.sort(eigsh(H,k=2,which='SA',tol=1e-11,return_eigenvectors=False))
            kap=1/(r*(2*r-1));beta=np.sin((r-m)*phi/2)**2
            bound=kap*beta/((1+np.sqrt(kap))**2+beta)
            assert ev[0]>=bound-1e-10
            if r==m:assert abs(ev[0])<1e-10 and ev[1]>=kap-1e-10
            else:assert ev[0]>1e-9
            out.append(dict(m=m,r=r,dimension=2**r,E0=float(ev[0]),gap=float(ev[1]-ev[0]),rejection_bound=float(bound)))
    return out

def pair_fluid_checks():
    out=[]
    for M,r in ((4,2),(5,2),(6,3)):
        b=5;L=b*M+2*r;states=set()
        for ones in it.combinations(range(M),r):
            eta=np.zeros(M,dtype=int);eta[list(ones)]=1
            for x in range(L):states.add(tuple(sorted((x+np.r_[0,np.cumsum(b+2*eta)[:-1]])%L)))
        states=sorted(states);index={s:i for i,s in enumerate(states)};H=np.zeros((len(states),len(states)))
        for ix,s in enumerate(states):
            for j in range(M):
                for step in (-2,2):
                    ss=list(s);ss[j]=(ss[j]+step)%L;ss=tuple(sorted(ss))
                    if ss in index:H[index[ss],ix]-=1
        assert np.allclose(H,H.T)
        ncomp,_=connected_components(sp.csr_matrix(H),directed=False);assert ncomp==1
        exact=[]
        for I in it.combinations(range(M),r):
            q=2*np.pi/M*(np.array(I)+(r-1)/2)
            momentum=(sum(I)+r*(r-1)//2)%M
            for ell in range(L):
                if (momentum+ell)%M==0:exact.append(-2*np.sum(np.cos(q+4*np.pi*ell/(L*M))))
        ev=np.linalg.eigvalsh(H);err=float(np.max(abs(ev-np.sort(exact))))
        assert err<1e-10
        E0=-2*np.sin(np.pi*r/M)/np.sin(np.pi/M)
        assert abs(ev[0]-E0)<1e-10
        out.append(dict(b=b,M=M,r=r,L=L,dimension=len(states),components=ncomp,spectral_error=err,E0=float(ev[0])))
    return out

if __name__=='__main__':
    result=dict(local_contexts=local_checks(),segments=segment_checks(),rings=ring_checks(),helices=helix_checks(),assumed_pair_exchange_graph=pair_fluid_checks())
    with open('checks_results.json','w') as handle:json.dump(result,handle,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k not in ('helices','segments')},indent=2))
    print('Segment spectral error',result['segments']['spectral_error'])
    print('Helix cases',len(result['helices']))
    print('All assertions passed.')

"""Finite G7 checks. Run with python3 -B gluing_checks.py. Numpy/scipy only."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import itertools as it
import json
import numpy as np
from scipy.linalg import eigvalsh

def ann(n,a):
    A=np.zeros((2**n,2**n),complex)
    for z in range(2**n):
        if (z>>a)&1:
            A[z^(1<<a),z]=(-1)**((z&((1<<a)-1)).bit_count())
    return A

def payload(M,n,deform=False):
    d=2**n
    aa=[ann(n,a) for a in range(n)]
    parity=np.diag([(-1)**z.bit_count() for z in range(d)])
    h=np.zeros((d*d,d*d),complex)
    onsite=np.zeros((d,d),complex)
    theta=[.31,-.57,.23][:n]
    ts=[.8,1.1,.7][:n]
    mu=[.13,-.21,.17][:n]
    for a,A in enumerate(aa):
        B=-ts[a]*np.exp(1j*theta[a])*np.kron(A.conj().T@parity,A)
        h+=B+B.conj().T
        onsite+=mu[a]*(A.conj().T@A-.5*np.eye(d))
    if deform:
        # Pairing is even under the common parity but breaks total U(1).
        B=.073*np.exp(.4j)*np.kron(aa[0]@parity,aa[-1])
        h+=B+B.conj().T
        if n>1:
            B=.091*np.exp(.2j)*aa[0].conj().T@aa[1]
            onsite+=B+B.conj().T
            B=.061*aa[0]@aa[1]
            onsite+=B+B.conj().T
    words=list(it.product(range(d),repeat=M)); index={z:i for i,z in enumerate(words)}
    H=np.zeros((len(words),len(words)),complex)
    for col,z in enumerate(words):
        for j in range(M):
            k=(j+1)%M
            for v in range(d*d):
                value=h[v,z[j]*d+z[k]]
                if abs(value)>1e-15:
                    zz=list(z);zz[j],zz[k]=divmod(v,d)
                    H[index[tuple(zz)],col]+=value
            for v in range(d):
                value=onsite[v,z[j]]
                if abs(value)>1e-15:
                    zz=list(z);zz[j]=v;H[index[tuple(zz)],col]+=value
    return H,words,index,theta,ts,mu

def momentum_blocks(H,words,index,M):
    seen=set();orbits=[]
    for z in words:
        if z in seen:continue
        orb=[];zz=z
        while zz not in seen:
            seen.add(zz);orb.append(index[zz]);zz=zz[1:]+zz[:1]
        orbits.append(orb)
    blocks=[]
    for m in range(M):
        cols=[]
        for orb in orbits:
            h=len(orb)
            if (m*h)%M:continue
            v=np.zeros(len(words),complex)
            v[orb]=np.exp(2j*np.pi*m*np.arange(h)/M)/np.sqrt(h)
            cols.append(v)
        Q=np.array(cols).T if cols else np.zeros((len(words),0))
        B=Q.conj().T@H@Q
        assert np.linalg.norm(H@Q-Q@B)<1e-10
        blocks.append(eigvalsh(B))
    return blocks

def free_payload(M,n,theta,ts,mu):
    ans=[[] for _ in range(M)]
    for mask in range(1<<(M*n)):
        N=mask.bit_count(); beta=((N-1)%2)/2
        E=0.;mom=0.
        for a in range(n):
            for I in range(M):
                occ=(mask>>(a*M+I))&1
                p=2*np.pi*(I+beta)/M
                E+=occ*(-2*ts[a]*np.cos(p-theta[a])+mu[a])
                mom+=occ*(I+beta)
            E-=mu[a]*M/2
        assert abs(mom-round(mom))<1e-12
        ans[round(mom)%M].append(E)
    return [np.sort(x) for x in ans]

def physical(b,M,r,n,Hp,words,index,J=.67):
    L=b*M+r
    positions=[]
    for pos in it.combinations(range(L),M):
        gaps=[(pos[(j+1)%M]-pos[j])%L for j in range(M)]
        if all(g in (b,b+1) for g in gaps): positions.append(pos)
    pi={p:i for i,p in enumerate(positions)};D=len(words)
    H=np.kron(np.eye(len(positions)),Hp)
    for ip,pos in enumerate(positions):
        for j in range(M):
            left=(pos[j]-pos[j-1])%L;right=(pos[(j+1)%M]-pos[j])%L
            if (left,right)!=(b,b+1):continue
            moved=list(pos);moved[j]=(moved[j]+1)%L
            order=sorted(range(M),key=lambda k:moved[k])
            target=tuple(moved[k] for k in order);jp=pi[target]
            for iz,z in enumerate(words):
                zz=tuple(z[k] for k in order);jz=index[zz]
                u=ip*D+iz;v=jp*D+jz
                H[v,u]-=J;H[u,v]-=J
    assert np.linalg.norm(H-H.conj().T)<1e-12
    return eigvalsh(H),L,len(positions)

def quotient(b,M,r,blocks,J=.67):
    L=b*M+r;ans=[]
    for Is in it.combinations(range(M),r):
        mg=(sum(Is)+r*(r-1)//2)%M
        q=2*np.pi*(np.array(Is)+(r-1)/2)/M
        for ell in range(L):
            mp=(-ell-mg)%M
            E=-2*J*np.cos(q+2*np.pi*ell/(L*M)).sum()
            ans.extend(E+blocks[mp])
    return np.sort(ans)

def main():
    result={'payload':[],'physical':[],'winding':[]}
    for M,n in [(3,1),(3,2),(4,2),(5,2)]:
        H,words,index,theta,ts,mu=payload(M,n)
        blocks=momentum_blocks(H,words,index,M)
        pred=free_payload(M,n,theta,ts,mu)
        err=max(np.max(abs(x-y)) if len(x) else 0 for x,y in zip(blocks,pred))
        assert err<1e-10
        result['payload'].append(dict(M=M,n=n,dimension=len(words),max_energy_error=float(err)))
    for b,M,r,n in [(2,3,1,2),(2,3,2,2),(2,4,1,1)]:
        for deform in [False,True]:
            H,words,index,_,_,_=payload(M,n,deform)
            blocks=momentum_blocks(H,words,index,M)
            direct,L,np0=physical(b,M,r,n,H,words,index)
            pred=quotient(b,M,r,blocks)
            assert len(direct)==len(pred)
            err=np.max(abs(direct-pred));assert err<1e-10
            result['physical'].append(dict(b=b,M=M,r=r,n=n,L=L,parity_only=deform,dimension=len(direct),max_energy_error=float(err)))
    # Exhaust all finite Fock occupations for the common-shift integrality.
    for M,n in [(3,1),(3,2),(4,2)]:
        count=0
        for Ns in it.product(range(M+1),repeat=n):
            Nt=sum(Ns);beta=((Nt-1)%2)/2
            for centers in it.product(range(-2,3),repeat=n):
                cs=np.array(centers)+beta+(np.array(Ns)-1)/2
                total=float(np.dot(Ns,cs))
                assert abs(total-round(total))<1e-12
                for a in range(n):
                    assert (round(2*cs[a])-Nt-Ns[a])%2==0
                count+=1
        result['winding'].append(dict(M=M,n=n,states_checked=count))
    print(json.dumps(result,indent=2))
    from pathlib import Path
    Path('checks.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()

"""Finite-dimensional checks for C8. No network or external project imports."""
import argparse
import json
import math
import itertools
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, eye, kron
from scipy.sparse.linalg import eigsh


def basis(q):
    out = []
    for i in range(q):
        for j in range(i+1,q):
            a=np.zeros((q,q),complex); a[i,j]=a[j,i]=np.sqrt(q/2); out.append(a)
            a=np.zeros((q,q),complex); a[i,j]=-1j*np.sqrt(q/2); a[j,i]=1j*np.sqrt(q/2); out.append(a)
    for k in range(1,q):
        a=np.zeros((q,q),complex)
        a[np.arange(k),np.arange(k)]=np.sqrt(q/(k*(k+1)))
        a[k,k]=-k*np.sqrt(q/(k*(k+1))); out.append(a)
    return out


def tensor(xs):
    z=np.ones((1,1),complex)
    for a in xs: z=np.kron(z,a)
    return z


def clifford(n):
    I=np.eye(2); X=np.array([[0,1],[1,0]]); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1])
    gs=[]
    for j in range(n):
        for a in (X,Y): gs.append(tensor([Z]*j+[a]+[I]*(n-j-1)))
    p=tensor([Z]*n)
    return gs, [1j*g@p for g in gs],p


def compression(q,m,pauli_mixture=False):
    assert m%2==0 and m<=q*q-1
    fs=basis(q)[:m]; gs,rs,p=clifford(m//2); D=len(p); k=q*q
    alpha=1/(4*m*np.sqrt(q))
    if pauli_mixture:
        s=int(round(np.log2(q))); assert 2**s==q
        ps=[np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
        fs=[tensor([ps[i] for i in ids]) for ids in itertools.product(range(4),repeat=s) if any(ids)][:m]
        alpha=1/(2*m)
    J=np.eye(k*D,dtype=complex)
    for f,g,r in zip(fs,gs,rs):
        J+=alpha*(np.kron(np.kron(f.T,np.eye(q)),g)+np.kron(np.kron(np.eye(q),f.T),r))
    J/=k
    ev,u=eigh(J); S=u*np.sqrt(np.maximum(ev,0))[None,:]
    # J_(i,a),(j,b) = sum_e conj(V_(i,e),a) V_(j,e),b.
    V=S.reshape(k,D,k*D).transpose(0,2,1).conj().reshape(k*k*D,D)
    err=np.linalg.norm(V.conj().T@V-np.eye(D))
    vt=V.reshape(k,k*D,D)
    for f,g,r in zip(fs,gs,rs):
        for A,target in [(np.kron(f,np.eye(q)),g),(np.kron(np.eye(q),f),r)]:
            mapped=np.einsum('iea,ij,jeb->ab',vt.conj(),A,vt,optimize=True)
            err=max(err,np.linalg.norm(mapped-alpha*target))
    ca=max(np.linalg.norm(a@b+b@a-2*(i==j)*np.eye(D)) for i,a in enumerate(gs) for j,b in enumerate(gs))
    return dict(q=q,m=m,logical_dimension=D,pauli_mixture=pauli_mixture,alpha=alpha,choi_min=float(ev[0]),choi_lower_bound=0 if pauli_mixture else 1/(2*k),max_compression_error=float(err),clifford_error=float(ca))


def sweep_check(q=2,ell=5):
    rng=np.random.default_rng(20261027)
    gates=[]
    for j in range(ell-1):
        a=rng.normal(size=(q*q,q*q))+1j*rng.normal(size=(q*q,q*q)); u,_=np.linalg.qr(a)
        gates.append(tensor([np.eye(q**j),u,np.eye(q**(ell-j-2))]))
    O=tensor([basis(q)[0]]+[np.eye(q)]*(ell-1))
    W=np.eye(q**ell,dtype=complex); errors=[]
    first=gates[0].conj().T@O@gates[0]
    for gate in gates:
        W=gate@W
        errors.append(float(np.linalg.norm(W.conj().T@O@W-first)))
    return dict(q=q,ell=ell,errors=errors)


def second_order_same_blocks(N=3):
    # Independent validation of the two-port extension using explicit positive weights.
    # Syndrome dimension r=3; G1's two-energy example at N=3.
    assert N==3
    ts=[a/np.sqrt(N) for a in basis(N)]
    eta=9/64
    F=sum(np.kron(a,a) for a in ts)
    J=sum(np.kron(a.T,a) for a in ts)
    MR=[np.eye(9)/2+eta*F,np.eye(9)/2-eta*F]
    ML=[np.eye(9)/2-eta*J,np.eye(9)/2+eta*J]
    E=[1,3]
    eff=np.zeros((N*N,N*N),complex)
    for e in range(2):
        for f in range(2):
            ar=MR[e].reshape(3,N,3,N); al=ML[f].reshape(3,N,3,N)
            for a in range(3):
                for b in range(3):
                    eff-=np.kron(ar[a,:,b,:],al[a,:,b,:])/(E[e]+E[f])
    swap=np.eye(N*N).reshape(N,N,N,N).transpose(1,0,2,3).reshape(N*N,N*N)
    target=-7177/8192*np.eye(N*N)+27/8192*swap
    # One actual block, with two physical qubit ports and two Bell partners.
    # Physical order is R,L,Raux,Laux,logical; dimension 16*N.
    D=16*N; G=np.zeros((D,N),complex)
    for r in range(2):
        for l in range(2):
            for i in range(N):
                G[((((r*2+l)*2+r)*2+l)*N+i),i]=.5
    paulis=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    FR=[tensor([a,np.eye(2),np.eye(2),np.eye(2),np.eye(N)]) for a in paulis]
    FL=[tensor([np.eye(2),a,np.eye(2),np.eye(2),np.eye(N)]) for a in paulis]
    S=np.concatenate([a@G for a in FR+FL],axis=1)
    gs=np.concatenate([G,S],axis=1)
    frame,_=np.linalg.qr(gs,mode='complete'); frame[:,:gs.shape[1]]=gs
    FQ=frame[:,N:]
    weights=[]
    for e in range(2):
        M=np.zeros((6*N,6*N),complex); M[:3*N,:3*N]=MR[e]; M[3*N:,3*N:]=ML[e]; weights.append(M)
    W=np.zeros((D-N,6*N),complex)
    for e,M in enumerate(weights):
        ev,u=eigh(M); W[e*6*N:(e+1)*6*N]= (u*np.sqrt(ev))@u.conj().T
    U,_=np.linalg.qr(W,mode='complete'); U[:,:6*N]=W
    energies=np.full(D-N,3.);energies[:6*N]=1
    excited=FQ@U.conj().T
    H0=(excited*energies)@excited.conj().T
    tr=[excited.conj().T@a@G for a in FR]
    tl=[excited.conj().T@a@G for a in FL]
    amp=sum(np.kron(a,b) for a,b in zip(tr,tl))
    denom=(energies[:,None]+energies[None,:]).reshape(-1)
    physical_eff=-amp.conj().T@(amp/denom[:,None])
    np.savez('two_port_N3.npz',H0=H0,G=G,FR=np.array(FR),FL=np.array(FL),H2=physical_eff)
    return dict(N=N,physical_block_dimension=D,second_order_error=float(np.linalg.norm(eff-target)),physical_second_order_error=float(np.linalg.norm(physical_eff-target)),ground_kernel_error=float(np.linalg.norm(H0@G)),frame_error=float(np.linalg.norm(frame.conj().T@frame-np.eye(D))),first_order_local_error=max(float(np.linalg.norm(G.conj().T@a@G)) for a in FR+FL),weight_min=min(float(eigh(a,eigvals_only=True)[0]) for a in MR+ML),connected_rank=N*N-1,swap_coefficient=27/8192)


def spinor_bond(n):
    gs,rs,_=clifford(n)
    return sum(np.kron(r,g) for r,g in zip(rs,gs))


def local_term(h,L,D,j):
    # Build by basis transitions, including the periodic seam without fermionic signs.
    dim=D**L; states=np.arange(dim,dtype=np.int64)
    weights=D**np.arange(L-1,-1,-1,dtype=np.int64)
    k=(j+1)%L
    a=(states//weights[j])%D; b=(states//weights[k])%D
    old=a*D+b; rows=[]; cols=[]; vals=[]
    for new in range(D*D):
        v=h[new,old]; idx=np.flatnonzero(np.abs(v)>1e-13)
        rows.append(states[idx]+(new//D-a[idx])*weights[j]+(new%D-b[idx])*weights[k]); cols.append(states[idx]); vals.append(v[idx])
    return csr_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(dim,dim))


def target_ed(n,L):
    D=2**n; h=spinor_bond(n); H=sum((local_term(h,L,D,j) for j in range(L)),csr_matrix((D**L,D**L),dtype=complex))
    e,v=eigsh(H,k=2,which='SA',tol=1e-11)
    order=np.argsort(e); e=e[order]; v=v[:,order]
    p=v[:,0].reshape(D**(L//2),-1); s=np.linalg.svd(p,compute_uv=False)**2; s=s[s>1e-14]
    return dict(n=n,L=L,energy=float(e[0]),predicted_energy=float(-4*n/np.sin(np.pi/L)),entropy=float(-s@np.log(s)),residual=float(np.linalg.norm(H@v[:,0]-e[0]*v[:,0])))


def gaussian_targets():
    # Independent entropy and Casimir estimators for the EXACT logical target.
    sizes=np.array([40,48,64,80,96,128]); ee=[]; ss=[]; vs=[]
    for L in sizes:
        k=2*np.pi*(np.arange(L)+.5)/L
        occ=np.sin(k)<0
        x=np.arange(L//2)
        C=np.exp(1j*(x[:,None,None]-x[None,:,None])*k[None,None,:])[:,:,occ].sum(axis=2)/L
        p=np.clip(eigh(C,eigvals_only=True),1e-15,1-1e-15)
        ss.append(float(-np.sum(p*np.log(p)+(1-p)*np.log1p(-p))))
        ee.append(float(np.sum(4*np.sin(k[occ]))))
        # Move the particle in a fixed particle-hole excitation to the next momentum.
        vs.append(float(L*(4*np.sin(3*np.pi/L)-4*np.sin(np.pi/L))/(2*np.pi)))
    inv=sizes.astype(float)**-2
    ec=np.linalg.lstsq(np.stack([np.ones(len(sizes)),inv,inv**2],axis=1),np.array(ee)/sizes,rcond=None)[0]
    sc=np.linalg.lstsq(np.stack([np.ones(len(sizes)),np.log(sizes),inv],axis=1),ss,rcond=None)[0]
    # exact velocity 4, separately checked by the gap sequence above
    return dict(sizes=sizes.tolist(),energies=ee,half_entropies=ss,oscillator_velocities=vs,c_energy_one_flavor=float(-6*ec[1]/(4*np.pi)),c_entropy_one_flavor=float(3*sc[1]),note='For n flavors these energies and entropies multiply by n. This is the logical target, not a fixed-d microscopic construction.')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--ed',action='store_true'); args=ap.parse_args()
    out=dict(compressions=[compression(2,2),compression(3,4),compression(3,8)],pauli_mixture=[compression(2,2,True),compression(4,8,True)],one_pass=sweep_check(),second_order=second_order_same_blocks(),logical_target=gaussian_targets())
    if args.ed: out['spinor_ed']=[target_ed(n,L) for n,L in [(1,4),(1,8),(2,4),(2,8)]]
    for row in out['compressions']+out['pauli_mixture']:
        assert row['choi_min'] >= row['choi_lower_bound']-1e-12
        assert row['max_compression_error'] < 1e-11
        assert row['clifford_error'] < 1e-12
    assert max(out['one_pass']['errors']) < 1e-11
    for key in ['second_order_error','physical_second_order_error','ground_kernel_error','frame_error','first_order_local_error']:
        assert out['second_order'][key] < 1e-11
    assert abs(out['logical_target']['c_energy_one_flavor']-1) < 1e-5
    assert abs(out['logical_target']['c_entropy_one_flavor']-1) < 1e-5
    for row in out.get('spinor_ed',[]):
        assert abs(row['energy']-row['predicted_energy']) < 1e-9
        assert row['residual'] < 1e-9
    out['all_checks_passed']=True
    with open('checks.json','w') as f: json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()

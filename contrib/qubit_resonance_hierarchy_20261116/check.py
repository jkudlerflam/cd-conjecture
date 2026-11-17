"""Exact-parent checks and full-sector higher-order checks with spectators."""
import json
import time
from itertools import combinations
from pathlib import Path
import numpy as np

from model import parent,qvalue,embed,hierarchy,maxabs,occupations


def ground_basis(L,m,r):
    q=qvalue(m);M=np.array([[-1.,1.],[-q,-1.]])/np.sqrt(1+q)
    S=np.diag([1.,r]);coeff=np.ones((1,1),complex)
    for j in range(L):
        V=S@np.linalg.matrix_power(M,j%m)
        d,n=coeff.shape;new=np.zeros((2*d,n+1),complex)
        for bit in (0,1):
            new[bit*d:(bit+1)*d,:n]+=V[bit,0]*coeff
            new[bit*d:(bit+1)*d,1:]+=V[bit,1]*coeff
        coeff=new
    U,sv,_=np.linalg.svd(coeff,full_matrices=False)
    assert min(sv)>max(sv)*1e-12
    return U,sv


def limiting_basis(L,m):
    cols=[];vac=np.zeros(1<<L);vac[0]=1;cols.append(vac)
    for residue in range(m):
        sites=list(range(residue,L,m))
        for N in range(1,len(sites)+1):
            v=np.zeros(1<<L)
            for occupied in combinations(sites,N):v[sum(1<<j for j in occupied)]=1
            cols.append(v/np.linalg.norm(v))
    return np.array(cols).T


def exact_parent_checks():
    rows=[];I=np.eye(2)
    for m in range(3,16):
        q=qvalue(m)
        for r in (.08,.23,.4):
            psi,h,G,T=parent(r,q)
            err=maxabs(h+q*q*r**4*np.eye(4)+np.kron(G,I)-np.kron(I,G)-np.outer(psi,psi.conj()))
            M=np.array([[-1.,1.],[-q,-1.]])/np.sqrt(1+q)
            power=np.linalg.matrix_power(M,m);phase_error=maxabs(power-np.trace(power)/2*np.eye(2))
            assert err<1e-12 and phase_error<1e-12
            rows.append({'m':m,'r':r,'bond_identity_error':err,'projective_period_error':phase_error})
    spectra=[]
    for m,L in ((3,6),(3,7),(4,8),(4,7),(5,10),(5,6),(6,6),(6,7)):
        psi,_,_,_=parent(.4,qvalue(m));H=embed(L,np.outer(psi,psi.conj()))
        ev=np.linalg.eigvalsh(H.toarray());count=int(sum(abs(ev)<1e-9));expected=L+1 if L%m==0 else 2
        assert count==expected and ev[count]>1e-8
        spectra.append({'m':m,'L':L,'ground_dimension':count,'first_positive_eigenvalue':float(ev[count])})
    limits=[];variational=[]
    for m in (3,4,5,6):
        L=2*m;G0=limiting_basis(L,m);previous=1.
        for r in (.1,.05,.025):
            G,_=ground_basis(L,m,r)
            sigma=np.linalg.svd(G0.conj().T@G,compute_uv=False)
            angle=float(np.sqrt(max(0.,1-min(sigma)**2)))
            assert angle<previous;previous=angle
            psi,_,_,_=parent(r,qvalue(m));H=embed(L,np.outer(psi,psi.conj()))
            err=maxabs(H@G);assert err<1e-10
            limits.append({'m':m,'L':L,'r':r,'subspace_distance':angle,'ground_residual':err})
        r=.23;q=qvalue(m);rho=r*np.sqrt(q)
        phi=np.array([1.,1j*rho])/np.sqrt(1+rho*rho);w=np.array([1j*rho,1.])/np.sqrt(1+rho*rho)
        k=2*np.pi/m+2*np.pi/L;trial=np.zeros(1<<L,complex)
        for j in range(L):
            v=np.ones(1,complex)
            for site in range(L):v=np.kron(w if site==j else phi,v)
            trial+=np.exp(1j*k*j)*v/np.sqrt(L)
        psi,_,_,_=parent(r,q);H=embed(L,np.outer(psi,psi.conj()));G,_=ground_basis(L,m,r)
        actual=float(np.vdot(trial,H@trial).real);pred=4*r*r*(1+q)*np.sin(np.pi/L)**2
        overlap=float(np.linalg.norm(G.conj().T@trial))
        assert abs(actual-pred)<1e-11 and overlap<1e-10
        variational.append({'m':m,'L':L,'r':r,'energy':actual,'predicted_bound':pred,'ground_overlap':overlap})
    return rows,spectra,limits,variational


def higher_case(m,L,periodic):
    P,K,stages=hierarchy(L,qvalue(m),periodic=periodic)
    M=K[2].toarray();pos={int(s):i for i,s in enumerate(P)};ns=occupations(L)[P]
    t=2048/(3*(123+55*np.sqrt(5))) if m==5 else 832/2187
    numbers=np.sum(ns,axis=1);adj=np.zeros_like(M)
    for col,state in enumerate(P):
        for j in range(L if periodic else L-m):
            k=(j+m)%L
            if ns[col,j]!=ns[col,k]:
                dest=state^(1<<j)^(1<<k)
                if dest in pos:adj[pos[dest],col]+=1
    degree=adj.sum(axis=0).real;potential=M+t*adj-np.diag(t*degree)
    off=potential.copy();np.fill_diagonal(off,0)
    residual=maxabs(off);pot=potential.diagonal().real
    mix=maxabs(M*(numbers[:,None]!=numbers[None,:]))
    assert residual<1e-9 and min(pot)>-1e-9 and mix<1e-9 and maxabs(K[1])<1e-9
    # Compare *every* single-residue sector with its coarse SWAP Laplacian.
    single_error=0.
    for label in range(m):
        indices=[i for i,s in enumerate(P) if all(j%m==label for j in np.flatnonzero(ns[i]))]
        single_error=max(single_error,maxabs(potential[np.ix_(indices,indices)]))
    assert single_error<1e-9
    eig=np.linalg.eigvalsh(M);kernel=int(sum(abs(eig)<1e-8))
    assert kernel==L+1
    # Check the exact sixth-constraint word kernel, not only its dimension.
    allowed=[]
    for state in range(1<<L):
        sites=[j for j in range(L) if (state>>j)&1]
        illegal=False
        for i,a in enumerate(sites):
            for b in sites[i+1:]:
                distances=(b-a,L-b+a) if periodic else (b-a,)
                if any(1<=d<=8 and d!=m for d in distances):illegal=True
        if not illegal:allowed.append(state)
    assert list(P)==allowed
    return {'m':m,'L':L,'periodic':periodic,'stages':stages,
            'dimension':len(P),'seventh_norm':maxabs(K[1]),'eighth_number_changing_norm':mix,
            'uniform_hop_graph_error':residual,'minimum_diagonal_potential':float(min(pot)),
            'single_residue_swap_error':single_error,'leading_kernel_dimension':kernel,
            'correct_constraint_words_verified':True,'hop_m':-t}


def main():
    start=time.perf_counter();a,b,c,d=exact_parent_checks()
    higher=[]
    for args in ((5,10,True),(5,15,True),(5,10,False),(5,13,False),(5,15,False),
                 (6,12,True),(6,12,False),(6,14,False),(6,16,False)):
        higher.append(higher_case(*args));print(higher[-1],flush=True)
    out={'numpy_version':np.__version__,'parent_identities':a,'ground_dimensions':b,
         'ground_space_limits':c,'quadratic_variational_bound':d,
         'higher_order_checks':higher,'runtime_seconds':time.perf_counter()-start}
    Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()

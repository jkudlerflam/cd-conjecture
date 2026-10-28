"""Finite-size diagnostic of a noncommuting, contact-coupled impurity.

H=Hgap+J*Delta/2*(1-sigma_x)+J*u*(A_0-alpha_L)*sigma_z.
The exact diagonalization uses the sigma_x eigenbasis. The perturbative
calculation retains J*Delta in its bath resolvent. It is not the d=24 model.
"""
from pathlib import Path
import argparse
import itertools
import json
import math
import os
import time
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla

HERE=Path(__file__).resolve().parent

def gap_chain(N,sites=(0,)):
    gaps=[]
    for cuts in itertools.combinations(range(2*N-1),N-1):
        markers=(-1,)+cuts+(2*N-1,)
        gaps.append(tuple(markers[j+1]-markers[j]-1 for j in range(N)))
    lookup={g:i for i,g in enumerate(gaps)}
    rows=[];cols=[]
    for col,g in enumerate(gaps):
        for j,value in enumerate(g):
            if value==0:continue
            for k in [(j-1)%N,(j+1)%N]:
                target=list(g);target[j]-=1;target[k]+=1
                rows.append(lookup[tuple(target)]);cols.append(col)
    dim=len(gaps)
    H=sp.csr_matrix((-np.ones(len(rows)),(rows,cols)),shape=(dim,dim))
    assert (H-H.T).nnz==0
    contacts=np.array([[float(g[site]==0) for site in sites] for g in gaps])
    A0=contacts[:,0]
    w,U=sla.eigsh(H,k=1,which='SA',v0=np.ones(dim),tol=2e-13)
    psi=U[:,0];E=float(w[0])
    if psi.sum()<0:psi=-psi
    alpha=float(np.dot(psi*psi,A0))
    exactE=-2/math.sin(math.pi/(2*N))
    exactAlpha=.5-1/(2*N*N*math.sin(math.pi/(2*N))**2)
    residual=float(np.linalg.norm(H@psi-E*psi))
    assert abs(E-exactE)<2e-11 and abs(alpha-exactAlpha)<2e-11
    assert residual<2e-10
    return H,psi,E,(A0 if len(sites)==1 else contacts),alpha,residual

def case(N,Delta=.2,u=1.):
    start=time.monotonic()
    H,psi,E0,A0,alpha,residual=gap_chain(N)
    dim=len(psi);F=A0-alpha;b=F*psi
    assert abs(np.dot(psi,b))<1e-12
    shifted=H-sp.eye(dim,format='csr')*E0
    spectral=None
    if N<=6:
        eigenvalues,eigenvectors=la.eigh(H.toarray())
        spectral=(eigenvalues-E0,abs(eigenvectors.T@b)**2)
    out=dict(N=N,physical_length=2*N,vacancies=N,dimension=dim,
             alpha_L=alpha,free_energy=E0,free_residual=residual,
             Delta=Delta,u=u,runs=[])
    for J in [.04,.02,.01]:
        mass=J*Delta
        resolvent=shifted+mass*sp.eye(dim,format='csr')
        response,info=sla.cg(resolvent,b,rtol=1e-12,atol=1e-14,maxiter=4000)
        response_residual=float(np.linalg.norm(resolvent@response-b))
        assert info==0 and response_residual<2e-10
        susceptibility=float(np.dot(b,response))
        spectral_value=None
        if spectral is not None:
            energies,weights=spectral
            spectral_value=float(np.sum(weights/(energies+mass)))
            assert abs(spectral_value-susceptibility)<2e-11
        off=sp.diags(J*u*F,format='csr')
        # sigma_x basis: the original sigma_z is now the off-diagonal Pauli.
        total=sp.bmat([[shifted,off],[off,shifted+mass*sp.eye(dim,format='csr')]],format='csr')
        initial=np.r_[psi,.07*psi];initial/=np.linalg.norm(initial)
        evals,states=sla.eigsh(total,k=2,which='SA',v0=initial,ncv=32,
                              tol=2e-12,maxiter=6000)
        state=states[:,0];energy=float(evals[0])
        eigen_residual=float(np.linalg.norm(total@state-energy*state))
        assert eigen_residual<2e-9
        gp,gm=state[:dim],state[dim:]
        second=-J*J*u*u*susceptibility
        record=dict(J=J,bare_gap=mass,energy_correction=energy,
                    correction_over_J2=energy/(J*J),second_order=second,
                    second_order_over_J2=-u*u*susceptibility,
                    relative_error_vs_second=(energy-second)/abs(second),
                    error_over_J4=(energy-second)/J**4,
                    dressed_gap=float(evals[1]-evals[0]),
                    sigma_x=float(np.dot(gp,gp)-np.dot(gm,gm)),
                    sigma_z=float(2*np.dot(gp,gm)),
                    F_sigma_z=float(2*np.dot(gp,F*gm)),
                    response_residual=response_residual,
                    eigenpair_residual=eigen_residual,
                    spectral_sum=spectral_value,
                    resolvent_sum=susceptibility)
        out['runs'].append(record)
        print(json.dumps(dict(stage='impurity',N=N,**record)),flush=True)
    out['seconds']=time.monotonic()-start
    return out

def pair_case(N,Delta=.2,u=1.):
    start=time.monotonic()
    H,psi,E0,contacts,alpha,residual=gap_chain(N,sites=(0,1))
    dim=len(psi);F=contacts-alpha
    shifted=H-sp.eye(dim,format='csr')*E0
    identity=sp.eye(dim,format='csr')
    out=dict(N=N,physical_length=2*N,dimension=dim,impurity_sites=[0,1],
             Delta=Delta,u=u,alpha_L=alpha,runs=[])
    projected=sla.LinearOperator((dim,dim),matvec=lambda x:shifted@x+psi*np.dot(psi,x),dtype=float)
    x,info=sla.cg(projected,F[:,1]*psi,rtol=1e-12,atol=1e-14,maxiter=4000)
    assert info==0
    cross_static=float(np.dot(F[:,0]*psi,x))
    out['static_cross_susceptibility']=cross_static
    out['predicted_asymptotic_error_over_J3']=-2*u**4*cross_static**2/Delta
    for J in [.04,.02,.01]:
        mass=J*Delta;A=shifted+mass*identity;responses=[]
        for j in range(2):
            b=F[:,j]*psi
            x,info=sla.cg(A,b,rtol=1e-12,atol=1e-14,maxiter=4000)
            assert info==0 and np.linalg.norm(A@x-b)<1e-10
            responses.append(x)
        susceptibility=np.array([[np.dot(F[:,i]*psi,responses[j]) for j in range(2)] for i in range(2)])
        blocks=[[None]*4 for _ in range(4)]
        for state in range(4):
            blocks[state][state]=shifted+mass*state.bit_count()*identity
            for j in range(2):
                target=state^(1<<j)
                blocks[target][state]=sp.diags(J*u*F[:,j],format='csr')
        total=sp.bmat(blocks,format='csr')
        initial=np.concatenate([psi,.07*psi,.053*psi,.01*psi])
        initial/=np.linalg.norm(initial)
        w,U=sla.eigsh(total,k=1,which='SA',v0=initial,ncv=36,tol=2e-12,maxiter=6000)
        energy=float(w[0]);residual=float(np.linalg.norm(total@U[:,0]-energy*U[:,0]))
        assert residual<2e-9
        second=-J*J*u*u*np.trace(susceptibility)
        record=dict(J=J,energy_correction=energy,correction_over_J2=energy/J**2,
                    second_order=second,second_order_over_J2=second/J**2,
                    relative_error_vs_second=(energy-second)/abs(second),
                    error_over_J3=(energy-second)/J**3,
                    susceptibility=susceptibility.tolist(),eigenpair_residual=residual)
        out['runs'].append(record)
        print(json.dumps(dict(stage='two_impurities',N=N,**record)),flush=True)
    out['seconds']=time.monotonic()-start
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-N',type=int,default=10)
    args=ap.parse_args();results=[]
    for N in [4,6,8,10]:
        if N>args.max_N:continue
        results.append(case(N))
        (HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    pairs=[]
    for N in [6,8,10]:
        if N>args.max_N:continue
        pairs.append(pair_case(N))
        (HERE/'pair_results.json').write_text(json.dumps(pairs,indent=2)+'\n')
    print(json.dumps(dict(stage='complete',seconds=sum(x['seconds'] for x in results+pairs))),flush=True)

if __name__=='__main__':main()

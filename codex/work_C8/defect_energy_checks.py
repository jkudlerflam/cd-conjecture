"""Checks for defect_energy_bounds.md; finite matrices and ring arithmetic."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,norm
from local_dressing_checks import embed


def case(M,seed=401):
    rng=np.random.default_rng(seed+M)
    d=5;b=4
    lengths=np.array([4,5,3,6,8]);energies=np.array([0.,0.,2.,3.,5.])
    Delta=2.;gamma=1.
    pair_length=(lengths[:,None]+lengths[None,:]).ravel()
    pg=((energies[:,None]+energies[None,:]).ravel()==0).astype(float)
    block=pg[:,None]*pg[None,:]+(1-pg[:,None])*(1-pg[None,:])
    dim=d**M;A=np.zeros((dim,dim),complex);W=np.zeros_like(A);F=np.zeros_like(A)
    Nt=np.zeros_like(A);Len=np.zeros_like(A);Tail=np.zeros_like(A)
    wi=np.zeros(M);fi=np.zeros(M)
    for j in range(M):
        A+=embed(np.diag(energies),j,1,M,d)
        Nt+=embed(np.diag(energies>0),j,1,M,d)
        Len+=embed(np.diag(lengths),j,1,M,d)
        Tail+=embed(np.diag(lengths>=7),j,1,M,d)
    for j in range(M-1):
        def random_term(scale,diagonal):
            z=rng.normal(size=(d*d,d*d))+1j*rng.normal(size=(d*d,d*d))
            z=(z+z.conj().T)/2
            z*=pair_length[:,None]==pair_length[None,:]
            if diagonal:z*=block
            return scale*z/norm(z,2)
        wd=random_term(.05,True);fl=random_term(.01,False)
        W+=embed(wd,j,2,M,d);F+=embed(fl,j,2,M,d)
        wi[j:j+2]+=.05;fi[j:j+2]+=.01
    w=wi.max();f=fi.max();gap=Delta-2*(w+f)
    ad=np.diag(A).real;nd=np.diag(Nt).real;ld=np.diag(Len).real;td=np.diag(Tail).real
    code=np.flatnonzero(ad==0)
    Eref=eigh((W+F)[np.ix_(code,code)],eigvals_only=True,subset_by_index=(0,0))[0]
    results=[]
    for constraint in (None,round(4.5*M)):
        inds=np.arange(dim) if constraint is None else np.flatnonzero(ld==constraint)
        pc=inds[ad[inds]==0]
        er=eigh((W+F)[np.ix_(pc,pc)],eigvals_only=True,subset_by_index=(0,0))[0]
        ev,U=eigh((A+W+F)[np.ix_(inds,inds)],subset_by_index=(0,0));psi=U[:,0];weights=abs(psi)**2
        E0=float(ev[0]);exc=max(0,float(er-Eref))
        eA=float(weights@ad[inds]);eN=float(weights@nd[inds]);tail=float(weights@td[inds])
        boundA=4*Delta*M*f*f/gap**2+2*Delta*exc/gap
        boundN=boundA/Delta;boundTail=boundA/(gamma*(7-b))
        assert E0>=Eref-M*f*f/gap-1e-12 and E0<=er+1e-12
        assert eA<=boundA and eN<=boundN and tail<=boundTail
        results.append(dict(total_length=constraint,dimension=len(inds),reference=float(er),
                            unrestricted_reference=float(Eref),ground=E0,
                            atomic_energy=eA,atomic_bound=boundA,defects=eN,defect_bound=boundN,
                            long_cells=tail,long_cell_bound=boundTail))
    return dict(cells=M,site_dimension=d,Delta=Delta,w=float(w),f=float(f),g=float(gap),sectors=results)


def ring_constraint_checks():
    results=[]
    for b in (3,9,39):
        worst=0.;count=0
        for M in range(3,101):
            for r in sorted({0,1,M//2,max(0,M//2-2),min(M,M//2+2),M-1,M}):
                L=b*M+r
                elen=-2*np.sin(np.pi*r/M)/np.sin(np.pi/M)
                minimum=-2*np.sin(np.pi*(M//2)/M)/np.sin(np.pi/M)
                for momentum in range(M):
                    ell=(-momentum)%M;phi=2*np.pi*ell/(L*M)
                    cost=elen*np.cos(phi)-minimum
                    bound=np.pi**2*(r-M/2)**2/(2*M)+2*np.pi**2/(b*b*M)
                    assert cost<=bound+1e-12
                    worst=max(worst,float(cost/bound));count+=1
        results.append(dict(b=b,cases=count,max_cost_over_bound=worst))
    return results


def main():
    out=dict(status='finite checks of energy, occupation and constraint bounds',
             atomic_cases=[case(M) for M in (2,3,4)],ring_constraint=ring_constraint_checks())
    Path('defect_energy_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()

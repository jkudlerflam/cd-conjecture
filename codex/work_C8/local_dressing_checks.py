"""Finite checks of local_dressing.md, with a two-dimensional local code."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, expm, norm, cho_factor, cho_solve


def embed(op, left, width, count, d):
    return np.kron(np.eye(d**left), np.kron(op, np.eye(d**(count-left-width))))


def fixture(count, high, seed=827):
    rng=np.random.default_rng(seed)
    d=4
    energies=np.array([0.,0.,1.,high])
    charges=np.array([0,1,0,1])
    pair_energy=(energies[:,None]+energies[None,:]).ravel()
    pair_charge=(charges[:,None]+charges[None,:]).ravel()
    pg=(pair_energy==0).astype(float); qg=1-pg
    rr=np.zeros(16);rr[qg>0]=1/pair_energy[qg>0]
    H0=np.zeros((d**count,d**count),complex)
    V=np.zeros_like(H0);S=np.zeros_like(H0);D=np.zeros_like(H0);O=np.zeros_like(H0)
    incident=np.zeros(count);generator_errors=[]
    for j in range(count):H0+=embed(np.diag(energies),j,1,count,d)
    for j in range(count-1):
        v=rng.normal(size=(16,16))+1j*rng.normal(size=(16,16))
        v=(v+v.conj().T)/2
        v*=pair_charge[:,None]==pair_charge[None,:]
        v/=norm(v,2)
        off=(qg[:,None]*pg[None,:]+pg[:,None]*qg[None,:])*v
        diag=v-off
        s=(rr[:,None]*pg[None,:]-pg[:,None]*rr[None,:])*v
        generator_errors.append(float(norm(s@np.diag(pair_energy)-np.diag(pair_energy)@s+off,2)))
        assert norm(s,2)<=1+1e-12
        assert norm(diag,2)<=1+1e-12 and norm(off,2)<=1+1e-12
        incident[j:j+2]+=1
        for total,part in ((V,v),(S,s),(D,diag),(O,off)):
            total+=embed(part,j,2,count,d)
    diag=np.diag(H0).real
    good=np.flatnonzero(diag==0);bad=np.flatnonzero(diag>0)
    R0=S@V-V@S-.5*(S@O-O@S)
    sw=-(V[np.ix_(good,bad)]/diag[bad])@V[np.ix_(bad,good)]
    second_error=float(norm(R0[np.ix_(good,good)]-sw,2))
    assert second_error<1e-11
    qtotal=sum(embed(np.diag(charges),j,1,count,d) for j in range(count))
    charge_error=float(norm(S@qtotal-qtotal@S,2))
    assert charge_error<1e-12
    J=float(incident.max())
    results=[]
    for eps in (.04,.02,.01):
        U=expm(eps*S)
        Hp=U@(H0+eps*V)@U.conj().T
        remainder=Hp-H0-eps*D
        remainder_norm=float(norm(remainder,2))
        bound=count*eps**2*12*J*J/(1-2*eps*J)**5
        assert remainder_norm<=bound
        expansion_error=float(norm(remainder/eps**2-R0,2))
        jr=12*J*J/(1-2*eps*J)**5
        gd=1-2*(eps*J+eps**2*jr)
        shift=-1-count*eps*J
        inverse=cho_solve(cho_factor(Hp-shift*np.eye(len(Hp)),lower=True),np.eye(len(Hp)))
        _,ground=eigh(inverse,subset_by_index=(len(Hp)-1,len(Hp)-1))
        e0=np.array([float(np.vdot(ground[:,0],Hp@ground[:,0]).real)])
        ground_residual=float(norm(Hp@ground[:,0]-e0[0]*ground[:,0]))
        er=float(eigh(Hp[np.ix_(good,good)],eigvals_only=True,subset_by_index=(0,0))[0])
        atomic=float(abs(ground[:,0])**2@diag)
        energy_bound=count*eps**4*jr**2/gd
        atomic_bound=4*count*eps**4*jr**2/gd**2
        assert gd>0 and -1e-12<=er-e0[0]<=energy_bound+1e-12 and atomic<=atomic_bound+1e-12
        results.append(dict(epsilon=eps,remainder_norm=remainder_norm,
                            rigorous_bound=bound,second_order_error=expansion_error,
                            dressed_energy_error=float(er-e0[0]),energy_error_bound=energy_bound,
                            dressed_atomic_energy=atomic,atomic_energy_bound=atomic_bound,
                            ground_residual=ground_residual))
    if high<1e4:
        assert results[-1]['second_order_error']<.6*results[-2]['second_order_error']
    # Slow time comparison, projected back to the code after dressing.
    observable=embed(np.diag([1,-1,0,0]),count//2,1,count,d)
    dynamics=[]
    for eps in (.04,.02,.01):
        U=expm(eps*S); Hp=U@(H0+eps*V)@U.conj().T
        ev,vec=eigh(Hp)
        T=.1;t=T/eps
        evolution=(vec*np.exp(-1j*t*ev))@vec.conj().T
        exact=evolution.conj().T@observable@evolution
        hcode=V[np.ix_(good,good)]
        ucode=expm(-1j*T*hcode)
        expected=ucode.conj().T@observable[np.ix_(good,good)]@ucode
        error=float(norm(exact[np.ix_(good,good)]-expected,2))
        jr=12*J*J/(1-2*eps*J)**5
        bound=6*eps*jr*T*np.exp(8*J*T)
        assert error<=bound
        dynamics.append(dict(epsilon=eps,slow_time=T,error=error,bound=float(bound)))
    return dict(cells=count,site_dimension=d,dimension=d**count,high_energy=high,
                interaction_norm=J,generator_error=max(generator_errors),
                second_order_identity_error=second_error,charge_error=charge_error,
                remainder=results,slow_dynamics=dynamics)


def main():
    cases=[fixture(n,h) for n,h in ((2,1),(3,1),(4,1),(3,17),(3,1e6))]
    out=dict(status='finite checks only; no CFT or all-order elimination',cases=cases)
    Path('local_dressing_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()

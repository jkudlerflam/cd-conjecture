"""Physical helix interpolation and local transport checks."""
import json
import numpy as np
from length_filter_checks import plus,tensors
from vacancy_parent import charge_hamiltonian


def helix_bond(phi):
    s=np.array([0,1,-np.exp(-1j*phi),0])/np.sqrt(2)
    return np.outer(s,s.conj())


def one_track(N,phi,pin_left):
    h=helix_bond(phi);a=np.eye(2)-np.outer(plus(0),plus(0).conj())
    H=np.zeros((2**N,)*2,complex)
    for j in range(N-1):H+=np.kron(np.kron(np.eye(2**j),h),np.eye(2**(N-j-2)))
    H+=np.kron(a,np.eye(2**(N-1))) if pin_left else np.kron(np.eye(2**(N-1)),a)
    G=tensors([plus(j*phi if pin_left else (j-N+1)*phi) for j in range(N)])
    e=np.linalg.eigvalsh(H)
    error=np.linalg.norm(H@G)
    assert error<1e-12 and e[1]>=1/(N*(2*N-1))-1e-12
    return float(error),float(e[1])


def G(N,phi):
    return tensors([np.kron(plus(j*phi),plus((j-N+1)*phi)) for j in range(N)])


def apply_at(v,op,site,total_sites,localdim):
    shape=(localdim,)*total_sites
    arr=np.moveaxis(v.reshape(shape),site,0)
    arr=(op@arr.reshape(localdim,-1)).reshape(arr.shape)
    return np.moveaxis(arr,0,site).reshape(-1)


def one_case(k,phi):
    b=k+1;U=1.;v=.1
    rs=[]
    for p in (1,2):
        holes,h,chi,e=charge_hamiltonian(k,p,U,v)
        rs.append(float(1-sum(c*c for s,c in zip(holes,chi) if 0 in s)))
    u1,u2=plus(-(b-1)*phi),plus(-b*phi)
    chi=np.sin(phi/2)**2
    ML=((np.eye(2)-np.outer(u2,u2.conj()))/rs[0]
        +(np.eye(2)-np.outer(u1,u1.conj()))/rs[1])/chi
    MR=ML.conj()
    expectations=[];error=0
    for r,N in zip(rs,(b,b+1)):
        l=plus(-(N-1)*phi);rr=plus((N-1)*phi)
        left=np.vdot(l,ML@l);right=np.vdot(rr,MR@rr)
        error=max(error,abs(r*left-1),abs(r*right-1))
        expectations.append([float(left.real),float(right.real)])
    assert np.linalg.eigvalsh(ML)[0]>-1e-12 and error<1e-11
    before=np.kron(G(b,phi),G(b+1,phi))
    after=np.kron(G(b+1,phi),G(b,phi))
    Ub=np.diag([1,np.exp(1j*b*phi)])
    op=np.exp(-1j*b*phi)*np.kron(Ub,Ub)
    transformed=apply_at(before,op,b,2*b+1,4)
    observed=np.vdot(after,transformed)
    beta=np.cos(phi/2)**(2*b)
    terr=abs(observed-beta)
    assert terr<1e-12
    parents=[one_track(N,phi,left) for N in (b,b+1) for left in (True,False)]
    # Check the actual compensated boundary matrices using one logical
    # qubit, whose original endpoint X/Y amplitudes are r/sqrt(k).
    X=np.array([[0,1],[1,0]],complex)
    projected=[r*np.vdot(plus(-(N-1)*phi),ML@plus(-(N-1)*phi))*X/np.sqrt(k)
               for r,N in zip(rs,(b,b+1))]
    bond0=np.kron(projected[0],projected[1])
    berr=np.linalg.norm(bond0-np.kron(X,X)/k)
    assert berr<1e-10
    return dict(k=k,b=b,phi=phi,r_values=rs,compensator_expectations=expectations,
                equalized_port_error=float(error),compensator_norm=float(np.linalg.norm(ML,2)),
                beta_helix=float(beta),local_conversion_error=float(terr),
                parent_max_residual=max(x[0] for x in parents),
                parent_min_gap=min(x[1] for x in parents),
                projected_two_cell_bond_error=float(berr),
                transport_state_dimension=len(before))


if __name__=='__main__':
    out=[one_case(3,np.pi/2),one_case(3,1/4)]
    with open('helix_equalizer_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))

"""Small independent algebra checks for C126; outputs only in this folder."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import expm,eigh
from checks import ground,bond,marginal,trial_parameters,ROOT

def ring_full(h,d,n,edges=None):
    dim=d**n;ids=np.arange(dim);powers=d**np.arange(n)
    H=np.zeros((dim,dim),complex)
    if edges is None:edges=[(x,(x+1)%n) for x in range(n)]
    for x,y in edges:
        a=ids//powers[x]%d;b=ids//powers[y]%d
        for ap in range(d):
            for bp in range(d):
                H[ids+(ap-a)*powers[x]+(bp-b)*powers[y],ids]+=h[ap*d+bp,a*d+b]
    return H

def cyclic(N):
    # The on-site spin shift on N qubits.
    codes=np.arange(2**N);dst=((codes<<1)&(2**N-1))+(codes>>(N-1))
    P=np.zeros((2**N,2**N),complex);P[dst,codes]=1
    return P

def half_cycle(N):
    P=cyclic(N);phi=2*np.pi*np.arange(-(N//2),N//2+1)/N
    out=np.zeros_like(P)
    # Principal logarithm, evaluated through spectral projectors of P.
    for angle in phi:
        proj=sum(np.exp(-1j*angle*r)*np.linalg.matrix_power(P,r) for r in range(N))/N
        out+=np.exp(.5j*angle)*proj
    return out

def run():
    out={'status':'numerical'}
    # Spin seam relabeling, including oriented complex bonds.
    rng=np.random.default_rng(126);X=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));h=(X+X.conj().T)/2
    errors=[]
    for L,N in [(2,2),(3,2),(2,3)]:
        ordinary=[];closing=[]
        for a in range(N):
            ordinary += [(a*L+x,a*L+x+1) for x in range(L-1)]
            closing += [(a*L+L-1,a*L)]
        codes=np.arange(2**(L*N));dst=codes.copy()
        for a in range(N):
            dst+=((codes>>(a*L))&1)*((1<<(((a+1)%N)*L))-(1<<(a*L)))
        U=np.zeros((len(codes),len(codes)));U[dst,codes]=1
        B=ring_full(h,2,L*N,closing)
        ht=ring_full(h,2,L*N,ordinary)+U@B@U.T
        hl=ring_full(h,2,L*N)
        errors.append(float(np.max(np.abs(ht-hl))))
    out['oriented_seam_matrix_errors']=errors
    # Full-ring ED cross-check of the balanced sectors on small rings.
    sectors=[]
    for d,n,delta in [(2,4,.5),(2,4,1.),(2,4,2.),(3,3,0.),(3,6,0.)]:
        e,*_=ground(d,n,delta);ef=float(eigh(ring_full(bond(d,delta),d,n),eigvals_only=True,subset_by_index=[0,0])[0])
        sectors.append(dict(d=d,L=n,delta=delta,error=abs(e-ef)))
    out['full_space_sector_checks']=sectors
    # Ground-state stabilizer and exact swap formula at a Heisenberg ring.
    e,c,v,res,dim=ground(2,6,1.);full=np.zeros(64);full[c]=v
    rho=marginal(full,2,6,2);sx=np.array([[0,1],[1,0]]);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1.]);I=np.eye(2)
    pauli=[I,sx,sy,sz];T=np.array([[np.trace(rho@np.kron(a,b)).real for b in pauli] for a in pauli])
    r=e/6;expected=np.diag([1.,4*r/3,4*r/3,4*r/3])
    out['Heisenberg_correlation_matrix_error']=float(np.max(abs(T-expected)))
    a,b,_,_,err=trial_parameters(full,2,6,bond(2,1.))
    out['Heisenberg_fractional_b_error']=abs(b-4*r*r/3)
    out['Heisenberg_two_copy_correlation_eigenvalues']=np.linalg.eigvalsh(np.kron(T,T)).tolist()
    # Uniform fractional swap has a strictly positive energy density.
    out['Heisenberg_uniform_halfway_excess_per_bond']=a+b
    # Independent full-vector expectation, with four nonuniform swap angles.
    e4,c4,v4,_,_=ground(2,4,1.);full4=np.zeros(16);full4[c4]=v4
    a4,b4,_,_,_=trial_parameters(full4,2,4,bond(2,1.))
    theta=np.array([.19,-.43,.72,1.13]);state=np.kron(full4,full4).astype(complex).reshape((2,)*8)
    for x,angle in enumerate(theta):
        state=np.cos(angle)*state+1j*np.sin(angle)*state.swapaxes(x,x+4)
    vec=state.reshape(-1);fullcost=float(np.vdot(vec,ring_full(bond(2,1.),2,8)@vec).real-2*e4)
    z=np.cos(2*theta);w=np.sin(2*theta);sgn=np.array([1.,1.,1.,-1.])
    functional=4*a4+np.sum(sgn*(-a4*z*np.roll(z,-1)+b4*w*np.roll(w,-1)))
    out['fractional_full_vector_error']=abs(fullcost-functional)
    out['fractional_full_vector_norm_error']=abs(np.vdot(vec,vec).real-1)
    # Prime-N cat state: compute local transformed spin expectations directly.
    cat=[]
    for N in [3,5,7]:
        U=half_cycle(N);codes=np.arange(2**N)
        magnet=np.array([1-2*((codes>>a)&1) for a in range(N)])
        moments=magnet@(abs(U)**2)
        energy=-float(np.sum(moments*moments)/(2**N*N))
        expect=-(N*N+2)/(3*N*N)
        cat.append(dict(N=N,energy_per_copy=energy,formula=expect,error=abs(energy-expect),unitarity_error=float(np.linalg.norm(U.conj().T@U-np.eye(2**N)))))
    out['cyclic_cat_checks']=cat
    # Half-filled XX, odd permutation phases, and the signed even cycle.
    free=[]
    for L in [4,6,10,30]:
        E=-1/np.sin(np.pi/L)
        for N in [2,3,5,9]:
            phi=2*np.pi*(np.arange(N)-(N-1)/2)/N
            exact=-1/np.sin(np.pi/(N*L))-N*E
            trial=E*(np.cos(phi/L).sum()-N)
            bound=np.pi**2*.5*(N-1/N)/(6*L)
            free.append(dict(L=L,N=N,exact=exact,trial=trial,error=abs(exact-trial),bound=bound))
    out['free_twist_checks']=free
    # Paired XY spin chain: build the signed fermion-copy trial as a full vector,
    # move its internal sign to the final closing bond, and test in the spin ring.
    L=4;n=2*L;Jx=1.;Jy=.3;field=.71
    hh=Jx*np.kron(sx,sx)/4+Jy*np.kron(sy,sy)/4+field*(np.kron(sz,I)+np.kron(I,sz))/4
    vals,vecs=eigh(ring_full(hh,2,L));psi=vecs[:,0];base=float(vals[0])
    codes=np.arange(2**n)
    par3=np.array([(-1)**int(k).bit_count() for k in range(2**L)])
    parity=int(np.rint(np.vdot(psi,par3*psi).real));boundary=-parity
    annih=np.array([[0.,1.],[0.,0.]])
    fs=[]
    for x in range(n):
        f=np.ones((1,1))
        for y in reversed(range(n)):
            f=np.kron(f,annih if y==x else (sz if y<x else I))
        fs.append(f)
    G=np.zeros((2**n,2**n),complex)
    for x in range(L):
        G+=(np.pi*x/(2*L))*(fs[x].T@fs[x+L]-fs[x+L].T@fs[x])
    trial=expm(G)@np.kron(psi,psi)
    sign=np.array([boundary**int(k>>L).bit_count() for k in codes])
    trial=sign*trial
    hlong=ring_full(hh,2,n)
    actual=float(np.vdot(trial,hlong@trial).real-2*base)
    hop=Jx*np.kron(sx,sx)/4+Jy*np.kron(sy,sy)/4
    ehop=float(np.vdot(psi,ring_full(hop,2,L)@psi).real/L)
    predicted=L*ehop*(2*np.cos(np.pi/(2*L))-2)
    defect=float(eigh(hlong,eigvals_only=True,subset_by_index=[0,0])[0]-2*base)
    out['paired_XY_signed_spin_trial']=dict(L=L,parity=parity,actual_cost=actual,predicted_cost=predicted,
        error=abs(actual-predicted),exact_defect=defect,
        total_parity=float(sum(abs(trial[k])**2*(-1)**int(k).bit_count() for k in codes)))
    # One-copy ordinary Ising twist and two-copy permutation twist.
    zz=np.kron(sz,sz);L=4
    ordinary=ring_full(-zz,2,L,[(i,i+1) for i in range(L-1)])+ring_full(zz,2,L,[(L-1,0)])
    perm=ring_full(-zz,2,2*L)
    out['Ising_ordered_checks']={'ordinary_Z2_twist_cost':float(np.linalg.eigvalsh(ordinary).min()+L),
        'permutation_twist_cost':float(np.linalg.eigvalsh(perm).min()+2*L),
        'product_cat_trial_cost':2.,'fixed_opposite_phase_sector_cost':4.}
    assert max(errors)<1e-12
    assert max(x['error'] for x in sectors)<1e-10
    assert out['Heisenberg_correlation_matrix_error']<1e-10
    assert out['fractional_full_vector_error']<1e-10
    assert max(x['error'] for x in cat)<1e-10
    assert max(x['error'] for x in free)<1e-10
    assert out['paired_XY_signed_spin_trial']['error']<1e-10
    assert actual>=defect-1e-10
    assert all(x['trial']<=x['bound']+1e-10 for x in free)
    (ROOT/'algebra_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['free_twist_checks','Heisenberg_two_copy_correlation_eigenvalues']},indent=2))
if __name__=='__main__':run()

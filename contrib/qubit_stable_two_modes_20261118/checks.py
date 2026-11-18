"""Independent bounded checks of the finite amplitude and diagnostic claims."""
from quartet import Sector, np, la, threadpool_limits
from diagnostics import calibrate
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import time


def occupation_check(L=12,N=4,momentum=0):
    sec=Sector(L,N,momentum)
    states=sorted(sec.lookup); lookup={s:i for i,s in enumerate(states)}
    embedding=np.zeros((len(states),len(sec.reps)),complex)
    for state,(col,shift) in sec.lookup.items():
        if col>=0:
            embedding[lookup[state],col]=np.exp(-2j*np.pi*momentum*shift/L)/np.sqrt(sec.periods[col])
    H=np.zeros((len(states),len(states)),complex)
    for col,state in enumerate(states):
        xs=[j for j in range(L) if (state>>j)&1]
        H[col,col]=sum(-2*((state>>((j+2)%L))&1)-(2/3)*((state>>((j+3)%L))&1)
                       for j in xs)
        for j in xs:
            for d,amp in ((1,1j),(-1,-1j),(2,-1),(-2,-1)):
                k=(j+d)%L
                if (state>>k)&1: continue
                new=state^(1<<j)^(1<<k)
                if new in lookup: H[lookup[new],col]+=amp
    direct=embedding.conj().T@H@embedding
    reduced=sec.hamiltonian().toarray()
    error=float(np.max(abs(direct-reduced)))
    orth=float(np.max(abs(embedding.conj().T@embedding-np.eye(embedding.shape[1]))))
    sigma=sec.sigma().toarray()
    comm=float(np.max(abs(sigma@reduced-reduced@sigma)))
    assert max(error,orth,comm)<1e-11
    return dict(L=L,N=N,momentum_index=momentum,dimension=len(states),
                reduced_dimension=len(sec.reps),reduction_error=error,
                isometry_error=orth,sigma_commutator_error=comm)


def quadratic_cancellation():
    # Four canonical fermion orbitals. Their total energy is zero and no
    # two-orbital denominator vanishes. Every pair coefficient is retained.
    energies=list(map(Fraction,[-3,-1]))+[Fraction(3,2),Fraction(5,2)]
    coeff={(i,j):Fraction(1+2*i+3*j,7) for i in range(4) for j in range(i+1,4)}
    def create(s,i,j):
        sign=1
        for k in (j,i):
            if (s>>k)&1: return None,0
            sign*=(-1)**((s&((1<<k)-1)).bit_count())
            s|=1<<k
        return s,sign
    terms=[]
    for (i,j),c in coeff.items():
        state,sign=create(0,i,j)
        den=energies[i]+energies[j]
        assert den!=0
        complement=[k for k in range(4) if not (state>>k)&1]
        out,sign2=create(state,*complement)
        assert out==15
        terms.append(-c*sign*coeff[tuple(complement)]*sign2/den)
    total=sum(terms)
    assert total==0 and sum(map(abs,terms))>0
    return dict(exact_sum=str(total),terms=list(map(str,terms)),
                absolute_sum=str(sum(map(abs,terms))))


def resolvent_check():
    mu=-9/4; z=(-mu-np.sqrt(mu*mu-4))/2
    exact=z**(np.arange(80)+1)/(2*(1-z/3))
    H=-2*mu*np.eye(80)-2*np.eye(80,k=1)-2*np.eye(80,k=-1)
    H[0,0]-=2/3
    rhs=np.eye(80)[:,0]; direct=la.solve(H,rhs,assume_a='pos')
    error=float(np.max(abs(direct[:24]-exact[:24])))
    assert error<1e-14
    # The exact inequality used to locate both decaying conjugate roots.
    q=Fraction(27,50)
    F=(1-q)*(3-q)*(1+q)**2/q**3
    assert F>16 and 3*q*q+4*q-3>0
    return dict(mu=mu,z=z,odd_halfline_error=error,
                binding_proof_rational_F=str(F),binding_proof_F_gt_16=True)


def pauli_check():
    X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.diag([1,-1]); I=np.eye(2); b=(X+1j*Y)/2; n=(I-Z)/2
    Delta=Fraction(3,7); V=Fraction(11,5); mu=Fraction(2,9); k=Fraction(4,3)
    Delta,V,mu,k=map(float,(Delta,V,mu,k))
    t=Delta**2/(2*V); w=1j*k*t; D=k*t/2; h=(V-mu)/2
    spin=Delta*(np.kron(X,X)-np.kron(Y,Y))/2+D*(np.kron(X,Y)-np.kron(Y,X))
    spin+=V*np.kron(Z,Z)/4-h*(np.kron(Z,I)+np.kron(I,Z))/2
    occ=V*np.kron(n,n)-mu*(np.kron(n,I)+np.kron(I,n))/2
    occ=occ.astype(complex)+w*np.kron(b,b.conj().T)+w.conjugate()*np.kron(b.conj().T,b)
    occ+=Delta*(np.kron(b,b)+np.kron(b.conj().T,b.conj().T))
    error=float(np.max(abs(spin-occ-(V/4-h)*np.eye(4))))
    assert error<1e-14
    return dict(two_site_Pauli_to_occupation_error=error)


def saved_result_checks():
    here=Path(__file__).resolve().parent
    q=json.loads((here/'quartet.json').read_text())
    micro=json.loads((here/'microscopic.json').read_text())
    exact=json.loads((here/'exact_dark.json').read_text())
    for r in q['rows']:
        assert r['intermediate_pole_distance']>.1
        assert r['spectral_sum_error']<1e-11 and r['resolvent_residual']<1e-11
        if r['connected']:
            assert r['connected']['identity_error']<1e-11
    assert exact['epsilon6_off_diagonal']=='6*sqrt(30)/5'
    assert exact['mu_eff_epsilon2']=='17/5'
    expected=12*np.sqrt(30)/5
    assert abs(micro['microscopic_dark_series']['predicted_gap_over_2_t_eps4']-expected)<1e-11
    assert micro['dark_fourth_order']['dark_pair_source_norm']<1e-11
    ratios=[r['ratio_to_prediction'] for r in micro['results'][0]['rows']]
    assert .93<ratios[0]<ratios[1]<ratios[2]<1
    assert ratios[-1]>.988
    darkvals=[r['coefficient_gap_over_2_t_eps4'] for r in micro['results'][2]['rows']]
    assert all(a<b for a,b in zip(darkvals,darkvals[1:])) and darkvals[-1]<expected
    # One Richardson removal of the epsilon^2 correction, not a phase fit.
    e1,e2=[r['epsilon'] for r in micro['results'][2]['rows'][-2:]]
    a1,a2=darkvals[-2:]
    extrap=(a2*e1**2-a1*e2**2)/(e1**2-e2**2)
    assert abs(extrap-expected)<.003
    return dict(generic_microscopic_ratios=ratios,
                dark_exact_gap_coefficient=expected,dark_last_two_Richardson=extrap,
                exact_to_numeric_series_error=abs(micro['microscopic_dark_series']['predicted_gap_over_2_t_eps4']-expected))


if __name__=='__main__':
    start=time.perf_counter(); here=Path(__file__).resolve().parent
    with threadpool_limits(limits=1):
        out=dict(occupation=[occupation_check(12,4,0),occupation_check(12,4,6)],
                 pauli=pauli_check(),free_fermion_quartet=quadratic_cancellation(),
                 exact_resolvent=resolvent_check(),saved_results=saved_result_checks(),
                 diagnostics=calibrate())
    sources=['codex/work_C179/PROOF.md','codex/work_C180/PROOF.md',
             'codex/work_C189/PROOF.md','codex/work_C191/PROOF.md','codex/work_C175/PROOF.md',
             'contrib/v4_locking_operators_20261118/AMPLITUDES.md',
             'contrib/v4_locking_operators_20261118/GAUSSIAN.md',
             'contrib/v4_locking_operators_20261118/LEDGER.md']
    root=here.parents[1]
    out['source_sha256']={s:hashlib.sha256((root/s).read_bytes()).hexdigest() for s in sources}
    out['seconds']=time.perf_counter()-start
    (here/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

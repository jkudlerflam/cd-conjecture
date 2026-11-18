"""Independent, bounded real-space validation of the sixth energy derivative.

No translation reduction is used in the Hamiltonian or resolvent. The
90-word spectral formula is compared to multivariate Rayleigh recursion.
"""
from medium import (Medium, GAMMA, REAL, IMAG, MU, np, sp, threadpool_limits)
from scipy.sparse.linalg import splu
from pathlib import Path
import itertools
import json
import time


def real_space_check(L=16):
    medium=Medium(L)
    states=sorted(s for sec in medium.bases.values() for s in sec.lookup)
    lookup={s:i for i,s in enumerate(states)}
    size=len(states); mask=(1<<L)-1
    def allowed(s):
        return not (s & (((s<<1)&mask)|(s>>(L-1))))
    def matrix(transitions):
        rows=[]; cols=[]; values=[]
        for col,s in enumerate(states):
            for out,v in transitions(s):
                if allowed(out):
                    rows.append(lookup[out]);cols.append(col);values.append(v)
        return sp.coo_matrix((values,(rows,cols)),shape=(size,size),dtype=complex).tocsc()
    def h_terms(s):
        diag=-MU*s.bit_count()
        for j in range(L):
            if not ((s>>j)&1): continue
            diag-=2*((s>>((j+2)%L))&1)+(2/3)*((s>>((j+3)%L))&1)
            for d,v in ((1,REAL+1j*IMAG),(-1,REAL-1j*IMAG),(2,-1),(-2,-1)):
                site=(j+d)%L
                if not ((s>>site)&1): yield s^(1<<j)^(1<<site),v
        yield s,diag
    def source_terms(s):
        for d,g in zip((2,3),GAMMA):
            for j in range(L):
                bits=(1<<j)|(1<<((j+d)%L))
                if not s&bits: yield s|bits,g
    def annihilation_terms(k):
        def transitions(s):
            for j in range(L):
                if (s>>j)&1:
                    yield s^(1<<j),np.exp(-2j*np.pi*k*j/L)/np.sqrt(L)
        return transitions
    H=matrix(h_terms)
    operators=[matrix(source_terms),matrix(annihilation_terms(-L//2+1)),
               matrix(annihilation_terms(-1))]
    # Only the normalized reference eigenvector uses the momentum calculation.
    # The full-space Hamiltonian checks both its energy and embedding phases.
    sec=medium.sector(medium.reference)
    ground=np.zeros(size,complex)
    for s,(i,shift) in sec.lookup.items():
        if i>=0:
            ground[lookup[s]]=medium.psi[i]*np.exp(-2j*np.pi*medium.m*shift/L)/np.sqrt(sec.periods[i])
    assert abs(np.vdot(ground,ground)-1)<1e-12
    ground_residual=float(np.linalg.norm(H@ground-medium.E0*ground))
    assert ground_residual<1e-10
    augmented=sp.bmat([[H-medium.E0*sp.eye(size),sp.csc_matrix(ground[:,None])],
                      [sp.csc_matrix(ground.conj()[None,:]),sp.csc_matrix((1,1))]],format='csc')
    factor=splu(augmented)
    residual=0.0
    def resolve(v):
        nonlocal residual
        rhs=np.r_[v,0]
        answer=factor.solve(rhs)
        residual=max(residual,float(np.linalg.norm(augmented@answer-rhs)/max(1,np.linalg.norm(rhs))))
        return answer[:-1]
    counts=sorted(itertools.product(range(3),repeat=3),key=lambda v:(sum(v),v))
    psi={(0,0,0):ground}; energies={}; target=(2,2,2)
    for nu in counts[1:]:
        rhs=np.zeros(size,complex)
        for axis in range(3):
            if nu[axis]:
                prev=list(nu);prev[axis]-=1
                rhs+=operators[axis]@psi[tuple(prev)]
        energies[nu]=np.vdot(ground,rhs)
        if nu==target: break
        # All proper energy counterterms vanish by charge and momentum.
        # Keep the general recursion anyway, so it does not assume that fact.
        for mu,e in energies.items():
            if mu!=nu and all(mu[j]<=nu[j] for j in range(3)):
                rest=tuple(nu[j]-mu[j] for j in range(3))
                rhs-=e*psi[rest]
        psi[nu]=-resolve(rhs)
    direct=energies[target]
    spectral=medium.four_leg()
    expected=complex(spectral['sixth_energy_coefficient_real'],spectral['sixth_energy_coefficient_imag'])
    proper=max(abs(e) for c,e in energies.items() if c!=target)
    error=abs(direct-expected)
    assert proper<1e-8
    assert error<1e-9*max(1,abs(expected))
    return dict(L=L,N=medium.N,full_dimension=size,ground_residual=ground_residual,
                largest_resolvent_relative_residual=residual,
                largest_proper_energy_counterterm=float(proper),
                real_space_sixth_coefficient=[float(direct.real),float(direct.imag)],
                spectral_sixth_coefficient=[float(expected.real),float(expected.imag)],
                absolute_difference=float(error),relative_difference=float(error/max(1e-30,abs(expected))),
                scope='L=16 validates the algorithm; its external momenta are outside the trial 1/8 patches used for L=20,24.')


def main():
    start=time.perf_counter()
    with threadpool_limits(limits=1): result=real_space_check()
    result['seconds']=time.perf_counter()-start
    Path(__file__).with_name('checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()

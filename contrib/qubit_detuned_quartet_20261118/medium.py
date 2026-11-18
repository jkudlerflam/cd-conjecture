"""Bounded finite-density source responses at the fixed p153c point.

The full leading interacting blockade Hamiltonian is retained. All number
and momentum sectors are searched to identify each finite reference ground
state; the physical chemical potential is never retuned. These finite
projected Legendre vertices are not thermodynamic continuum C4 estimates.
"""
from __future__ import annotations

import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'qubit_stable_two_modes_20261118'))
from quartet import Sector
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from threadpoolctl import threadpool_limits
import itertools
import json
import time
import argparse

PHI=1.53
EPS=1/40
T=1/80
MU=-2.1426
REAL=float(np.cos(PHI))
IMAG=float(np.sin(PHI))
GAMMA=np.array([-2*EPS*REAL,4*EPS/3])


def legal_bases(L):
    """Generate the Fibonacci-sized constrained space, without scanning 2^L."""
    groups=[[] for _ in range(L//2+1)]
    def walk(j,state,previous):
        if j==L:
            if not (previous and state&1): groups[state.bit_count()].append(state)
            return
        walk(j+1,state,False)
        if not previous: walk(j+1,state|(1<<j),True)
    walk(0,0,False)
    result={}
    for N,states in enumerate(groups):
        sec=object.__new__(Sector)
        sec.L,sec.N,sec.momentum,sec.blockade=L,N,0,True
        sec.mask=(1<<L)-1; sec.reps=[]; sec.periods=[]; sec.lookup={}
        for state in sorted(states):
            if state in sec.lookup: continue
            orbit=[]; s=state
            while s not in orbit:
                orbit.append(s); s=sec.translate(s)
            idx=len(sec.reps)
            sec.reps.append(state); sec.periods.append(len(orbit))
            for shift,s in enumerate(orbit): sec.lookup[s]=(idx,shift)
        sec.periods=np.array(sec.periods)
        result[N]=sec
    return result


def momentum_sector(base,m):
    sec=object.__new__(Sector)
    sec.L,sec.N,sec.momentum,sec.blockade=base.L,base.N,m,True
    sec.mask=base.mask
    indices=np.flatnonzero((m*base.periods)%base.L==0)
    indexmap={int(old):new for new,old in enumerate(indices)}
    sec.reps=[base.reps[i] for i in indices]
    sec.periods=base.periods[indices]
    sec.lookup={s:(indexmap.get(i,-1),shift) for s,(i,shift) in base.lookup.items()}
    return sec


class Medium:
    def __init__(self,L):
        self.L=L; self.bases=legal_bases(L)
        self.sectors={}; self.hs={}; self.spectra={}; self.pairs={}; self.singles={}
        records=[]; best=None; largest=0
        for N,base in self.bases.items():
            for m in range(L):
                sec=momentum_sector(base,m)
                if not sec.reps: continue
                H=sec.hamiltonian(IMAG,REAL)-MU*N*sp.eye(len(sec.reps))
                largest=max(largest,len(sec.reps))
                if len(sec.reps)==1:
                    es=np.array([H[0,0].real]); vs=np.ones((1,1),complex)
                elif len(sec.reps)<180:
                    es,vs=la.eigh(H.toarray(),subset_by_index=(0,1))
                else:
                    es,vs=eigsh(H,k=2,which='SA',tol=2e-12,
                                v0=np.ones(H.shape[0],complex))
                    order=np.argsort(es); es=es[order]; vs=vs[:,order]
                records.extend((float(e),N,m) for e in es)
                if best is None or es[0]<best[0]:
                    best=(float(es[0]),N,m,vs[:,0],sec,H)
        records.sort()
        self.E0,self.N,self.m,self.psi,sec,H=best
        self.reference=(self.N,self.m)
        self.sectors[self.reference]=sec; self.hs[self.reference]=H
        self.gap=records[1][0]-records[0][0]
        assert self.gap>1e-8, 'This calculation requires a unique finite ground state.'
        self.reference_info=dict(L=L,N=self.N,density=self.N/L,
            translation_momentum_index=self.m,physical_momentum_over_pi=-2*self.m/L,
            grand_energy_over_t=self.E0,global_gap_over_t=self.gap,
            first_levels=[dict(energy=e,N=n,momentum_index=m) for e,n,m in records[:8]],
            max_searched_sector=largest,total_constrained_dimension=sum(len(s.lookup) for s in self.bases.values()))

    def sector(self,key):
        N,m=key
        if N<0 or N>self.L//2: return None
        if key not in self.sectors:
            self.sectors[key]=momentum_sector(self.bases[N],m)
        return self.sectors[key]

    def spectrum(self,key):
        if key not in self.spectra:
            sec=self.sector(key)
            if key not in self.hs:
                self.hs[key]=sec.hamiltonian(IMAG,REAL)-MU*sec.N*sp.eye(len(sec.reps))
            es,vs=la.eigh(self.hs[key].toarray())
            delta=es-self.E0
            assert min(delta,default=1)>-1e-8
            self.spectra[key]=(delta,vs)
        return self.spectra[key]

    def resolve(self,key,vec,cutoff=0.0,power=1):
        if key is None: return vec
        delta,U=self.spectrum(key)
        keep=(delta>1e-9)&(delta>=cutoff)
        inverse=np.zeros(len(delta)); inverse[keep]=delta[keep]**(-power)
        return U@(inverse*(U.conj().T@vec))

    def creation(self,key,d,qindex):
        target=(key[0]+2,(key[1]-qindex)%self.L)
        src=self.sector(key); dst=self.sector(target)
        if src is None or dst is None: return target,None
        cache=(key,d,qindex)
        if cache not in self.pairs:
            def transitions(state):
                for j in range(self.L):
                    bits=(1<<j)|(1<<((j+d)%self.L))
                    if state&bits: continue
                    out=state|bits
                    if dst.allowed(out):
                        phase=np.exp(2j*np.pi*qindex*j/self.L)
                        yield out,phase
            self.pairs[cache]=src.matrix(transitions,dst)
        return target,self.pairs[cache]

    def act(self,op,key,vec,weights):
        """B/b create/annihilate a uniform source; P/p annihilate/create pi pair."""
        if key is None: return None,np.zeros(0,complex)
        source=op in ('B','b'); creation=op in ('B','p')
        q=0 if source else self.L//2
        target=(key[0]+(2 if creation else -2),(key[1]+q)%self.L)
        dst=self.sector(target)
        if dst is None or not dst.reps: return None,np.zeros(0,complex)
        out=np.zeros(len(dst.reps),complex)
        coeffs=zip((2,3),weights) if source else [(2,1/np.sqrt(self.L))]
        for d,coeff in coeffs:
            if creation:
                _,mat=self.creation(key,d,q)
                out+=coeff*(mat@vec)
            else:
                _,mat=self.creation(target,d,q)
                out+=np.conj(coeff)*(mat.conj().T@vec)
        return target,out

    def single_creation(self,key,k):
        target=(key[0]+1,(key[1]-k)%self.L)
        src=self.sector(key); dst=self.sector(target)
        if src is None or dst is None: return target,None
        cache=(key,k)
        if cache not in self.singles:
            def transitions(state):
                for j in range(self.L):
                    if (state>>j)&1: continue
                    out=state|(1<<j)
                    if dst.allowed(out):
                        yield out,np.exp(2j*np.pi*k*j/self.L)/np.sqrt(self.L)
            self.singles[cache]=src.matrix(transitions,dst)
        return target,self.singles[cache]

    def single_act(self,k,creation,key,vec):
        if key is None: return None,np.zeros(0,complex)
        target=(key[0]+(1 if creation else -1),
                (key[1]+(-k if creation else k))%self.L)
        dst=self.sector(target)
        if dst is None or not dst.reps: return None,np.zeros(0,complex)
        if creation:
            _,mat=self.single_creation(key,k)
            out=mat@vec
        else:
            _,mat=self.single_creation(target,k)
            out=mat.conj().T@vec
        return target,out

    def four_leg(self,cutoff=0.0,weights=GAMMA):
        """Static four-single-leg vertex of the two-mode Legendre transform.

        External physical momenta kA=-pi+2pi/L and kB=-2pi/L have sum -pi.
        No pair of external annihilation momenta sums to zero; no proper
        subset of BBaabb is charge- and momentum-neutral. Thus every Wick,
        folded-energy, and anomalous-leg subtraction vanishes exactly in
        this particular projection. The sixth derivative contains 90 words.
        """
        modes={'a':(-self.L//2+1)%self.L,'b':(-1)%self.L}
        for k1,k2 in itertools.product(modes.values(),repeat=2):
            assert (k1+k2)%self.L!=0
        def apply(op,key,vec):
            if op=='B': return self.act('B',key,vec,weights)
            if op=='C': return self.act('b',key,vec,weights)
            return self.single_act(modes[op.lower()],op.isupper(),key,vec)
        cache={(): (self.reference,self.psi)}
        def resolved_prefix(word):
            word=tuple(word)
            if word not in cache:
                key,v=resolved_prefix(word[:-1]); key,v=apply(word[-1],key,v)
                if key is not None: v=self.resolve(key,v,cutoff)
                cache[word]=(key,v)
            return cache[word]
        def three(word):
            key,v=resolved_prefix(word[:2])
            return apply(word[2],key,v)
        # Capital B denotes the source; use an explicit adjoint application
        # for the valley-b leg to avoid any label collision.
        def adjoint_three(word):
            key,v=self.reference,self.psi
            for i,op in enumerate(reversed(word)):
                if op=='B': key,v=self.act('b',key,v,weights)
                else: key,v=self.single_act(modes[op],True,key,v)
                if key is not None and i<2: v=self.resolve(key,v,cutoff)
            return key,v
        words=sorted(set(itertools.permutations('BBaabb')))
        terms=[]; c=0j; bound=0.0
        for word in words:
            kr,z=three(word[:3]); kl,y=adjoint_three(word[3:])
            if kr is None or kl is None:
                val=0j; cap=0.0
            else:
                assert kr==kl
                rz=self.resolve(kr,z,cutoff); ry=self.resolve(kl,y,cutoff)
                val=-np.vdot(y,rz)
                cap=float(np.sqrt(max(0,np.vdot(z,rz).real)*max(0,np.vdot(y,ry).real)))
                assert abs(val)<=cap+1e-8
            c+=val; bound+=cap
            terms.append(dict(word=''.join(word),real=float(val.real),imag=float(val.imag),bound=cap))
        chis=[]
        for k in modes.values():
            chi=0.0
            for dagger in (False,True):
                key,v=self.single_act(k,dagger,self.reference,self.psi)
                if key is not None: chi+=float(np.vdot(v,self.resolve(key,v,cutoff)).real)
            assert chi>0; chis.append(chi)
        amp=self.L*c/(chis[0]**2*chis[1]**2)
        amp_bound=self.L*bound/(chis[0]**2*chis[1]**2)
        return dict(cutoff_over_t=cutoff,number_of_words=len(words),
                    external_k_over_pi=[-1+2/self.L,-2/self.L],
                    single_leg_susceptibilities_times_t=chis,
                    sixth_energy_coefficient_real=float(c.real),
                    sixth_energy_coefficient_imag=float(c.imag),
                    spectral_Cauchy_bound=bound,
                    g4_L_over_t=[float(amp.real),float(amp.imag)],
                    g4_L_Cauchy_bound_over_t=float(amp_bound),
                    source_weights=np.asarray(weights).tolist(),
                    subtractions=dict(Wick=0,folded_energy=0,anomalous_external_legs=0),
                    terms=terms,
                    scope='Finite static 1PI coefficient in the Legendre transform of the two selected complex modes; g4_L=L times the collective coefficient. This normalization does not establish locality or thermodynamic matching.')

    def first_two(self,ops,weights,cutoff):
        key,vec=self.act(ops[0],self.reference,self.psi,weights)
        if key is None: return key,vec
        vec=self.resolve(key,vec,cutoff)
        return self.act(ops[1],key,vec,weights)

    def coefficient(self,weights=GAMMA,cutoff=0.0):
        weights=np.asarray(weights)
        adj={'B':'b','b':'B','P':'p','p':'P'}
        terms=[]; bounds=[]
        words=sorted(set(itertools.permutations('BBPP')))
        for word in words:
            kr,z=self.first_two(word[:2],weights,cutoff)
            kl,y=self.first_two([adj[word[3]],adj[word[2]]],weights,cutoff)
            if kr is None or kl is None:
                term=0j; bound=0.0
            else:
                assert kr==kl
                Rz=self.resolve(kr,z,cutoff); Ry=self.resolve(kl,y,cutoff)
                term=-np.vdot(y,Rz)
                bound=float(np.sqrt(max(0,np.vdot(z,Rz).real)*max(0,np.vdot(y,Ry).real)))
                assert abs(term)<=bound+1e-9
            terms.append(dict(word=''.join(word),real=float(term.real),imag=float(term.imag),bound=bound))
            bounds.append(bound)
        # Ground-state folded/disconnected terms vanish: a uniform source and
        # a momentum-pi probe cannot return the reference after two insertions.
        c=sum(complex(x['real'],x['imag']) for x in terms)
        chi=0.0
        spectral=[]
        for op in ('P','p'):
            key,v=self.act(op,self.reference,self.psi,weights)
            if key is None: continue
            rv=self.resolve(key,v,cutoff)
            chi+=float(np.vdot(v,rv).real)
            delta,U=self.spectrum(key)
            sw=abs(U.conj().T@v)**2
            spectral.append(dict(charge=2 if op=='p' else -2,
                                  total_weight=float(sum(sw)),
                                  retained_weight=float(sum(sw[(delta>1e-9)&(delta>=cutoff)])),
                                  smallest_weighted_gap=float(min(delta[sw>1e-10],default=np.inf))))
        assert chi>0
        return dict(cutoff_over_t=cutoff,source_weights=weights.tolist(),
                    coefficient_real=float(c.real),coefficient_imag=float(c.imag),
                    coefficient_abs=float(abs(c)),source_projected_bound=sum(bounds),
                    pair_susceptibility_times_t=chi,
                    composite_inverse_response_abs_over_t=float(abs(c)/chi**2),
                    composite_inverse_response_bound_over_t=float(sum(bounds)/chi**2),
                    disconnected_counterterm=0,terms=terms,pair_spectral_weights=spectral)

    def run(self):
        rows=[self.coefficient(cutoff=x) for x in (0.0,.25,.5)]
        channels=[self.coefficient(w) for w in ([1,0],[0,1],[1,1])]
        c=[complex(x['coefficient_real'],x['coefficient_imag']) for x in channels]
        cross=(c[2]-c[0]-c[1])/2
        reconstructed=GAMMA[0]**2*c[0]+2*GAMMA[0]*GAMMA[1]*cross+GAMMA[1]**2*c[1]
        measured=complex(rows[0]['coefficient_real'],rows[0]['coefficient_imag'])
        assert abs(reconstructed-measured)<1e-10
        return dict(reference=self.reference_info,vertices=rows,
                    source_channel_coefficients=dict(B2B2=[c[0].real,c[0].imag],
                         B2B3=[cross.real,cross.imag],B3B3=[c[1].real,c[1].imag]),
                    interference_reconstruction_error=float(abs(reconstructed-measured)))

    def direct_probe_check(self,taus):
        """Independent finite-field check of the mixed fourth derivative.

        The source and the pi probe close on two momentum sectors and one
        parity. No field, density, or physical coupling is retuned with L.
        """
        keys=[(N,m) for N in self.bases if N%2==self.N%2
              for m in (self.m,(self.m+self.L//2)%self.L)
              if self.sector((N,m)).reps]
        sizes=[len(self.sector(k).reps) for k in keys]
        ends=np.cumsum([0]+sizes); slices={k:slice(ends[i],ends[i+1]) for i,k in enumerate(keys)}
        dim=int(ends[-1]); H=np.zeros((dim,dim),complex)
        X=np.zeros_like(H); P=np.zeros_like(H)
        for key in keys:
            cut=slices[key]
            sec=self.sector(key)
            h=sec.hamiltonian(IMAG,REAL)-MU*sec.N*sp.eye(len(sec.reps))
            H[cut,cut]=h.toarray()
            for d,g in zip((2,3),GAMMA):
                target,mat=self.creation(key,d,0)
                if target in slices: X[slices[target],cut]+=g*mat.toarray()
            target,mat=self.creation(key,2,self.L//2)
            if target in slices: P[slices[target],cut]+=mat.toarray()/np.sqrt(self.L)
        X+=X.conj().T
        target_c=self.coefficient()['coefficient_real']
        rows=[]
        for tau in taus:
            energies=[]
            for phase in (1,1j):
                h=H+tau*X+tau*(phase*P+np.conj(phase)*P.conj().T)
                energies.append(float(la.eigvalsh(h,subset_by_index=(0,0))[0]))
            estimate=(energies[0]-energies[1])/(4*tau**4)
            rows.append(dict(tau=tau,mixed_derivative_estimate=estimate,
                             exact_resolvent_coefficient=target_c,
                             relative_error=(estimate/target_c-1)))
        return dict(dimension=dim,rows=rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['all','pair','four-leg'],default='all')
    args=parser.parse_args()
    start=time.perf_counter()
    with threadpool_limits(limits=1):
        results=[]; four=[]
        for L in (16,20,24):
            m=Medium(L)
            if args.mode in ('all','pair'):
                out=m.run(); results.append(out)
                if L==16: out['direct_probe_check']=m.direct_probe_check([1/8,1/16,1/32])
            if args.mode in ('all','four-leg') and L in (20,24):
                rows=[m.four_leg(cutoff=x) for x in (0.0,.25,.5)]
                channels=[m.four_leg(weights=w) for w in ([1,0],[0,1],[1,1])]
                cc=[complex(x['sixth_energy_coefficient_real'],x['sixth_energy_coefficient_imag']) for x in channels]
                cross=(cc[2]-cc[0]-cc[1])/2
                recon=GAMMA[0]**2*cc[0]+2*GAMMA[0]*GAMMA[1]*cross+GAMMA[1]**2*cc[1]
                actual=complex(rows[0]['sixth_energy_coefficient_real'],rows[0]['sixth_energy_coefficient_imag'])
                assert abs(recon-actual)<1e-9*max(1,abs(actual))
                four.append(dict(reference=m.reference_info,vertices=rows,
                    source_channel_coefficients=dict(B2B2=[cc[0].real,cc[0].imag],
                        B2B3=[cross.real,cross.imag],B3B3=[cc[1].real,cc[1].imag]),
                    B3_only_sixth_coefficient=[float((GAMMA[1]**2*cc[1]).real),float((GAMMA[1]**2*cc[1]).imag)],
                    interference_reconstruction_error=float(abs(recon-actual))))
            print('completed',L,'N',m.N,'gap',m.gap,flush=True)
    result=dict(parameters=dict(phi=PHI,Delta=1,V=40,t=T,epsilon=EPS,
                    real_hop_over_t=REAL,imaginary_hop_over_t=IMAG,
                    chemical_potential_over_t=MU,pair_sources_over_t=GAMMA.tolist()),
                results=results,seconds=time.perf_counter()-start,
                scope='Finite interacting composite-pair diagnostic, not single-field quartet 1PI or continuum C4.')
    if results: Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    if four:
        four_result=dict(parameters=result['parameters'],results=four,
                         seconds=time.perf_counter()-start,
                         scope='Finite interacting four-single-leg projected Legendre coefficient; not continuum C4.')
        Path(__file__).with_name('four_leg.json').write_text(json.dumps(four_result,indent=2)+'\n')
    print('seconds',time.perf_counter()-start)


if __name__=='__main__': main()

"""Exact spectra and physical-marker entropies for the projected length fluid.

This is a constrained marker model, not the finite-coupling microscopic
selector/wrapper Hamiltonian. Central charge estimates refer only to it.
"""
import argparse
import itertools
import json
import numpy as np
from scipy import sparse as sp
from scipy.linalg import svdvals


def marker_states(b,M,r):
    """Every physical marker pattern once, with its gap variables."""
    for longs in itertools.combinations(range(M),r):
        gaps=np.full(M,b,dtype=int);gaps[list(longs)]+=1
        relative=np.r_[0,np.cumsum(gaps)[:-1]]
        for first in range(gaps[-1]):
            yield tuple(map(int,relative+first)),longs


def physical_graph(b,M,r):
    L=b*M+r
    states=sorted(s for s,_ in marker_states(b,M,r))
    index={s:i for i,s in enumerate(states)}
    assert len(index)==len(states)
    H=np.zeros((len(states),)*2)
    for i,s in enumerate(states):
        for j,x in enumerate(s):
            if (x-s[j-1])%L==b and (s[(j+1)%M]-x)%L==b+1:
                target=list(s);target[j]=(x+1)%L;target=tuple(sorted(target))
                out=index[target];H[out,i]-=1;H[i,out]-=1
    return H


def exact_levels(b,M,r):
    L=b*M+r
    levels=[]
    for occupied in itertools.combinations(range(M),r):
        m=(sum(occupied)+r*(r-1)//2)%M
        ks=2*np.pi*(np.asarray(occupied)+(r-1)/2)/M
        for ell in range((-m)%M,L,M):
            P=2*np.pi*ell/L
            levels.append(float(-2*np.cos(ks+P/M).sum()))
    return np.sort(levels)


def spectral_checks():
    out=[]
    for b,M,r in [(2,5,2),(2,5,3),(3,5,2),(4,4,2),(4,6,3),(4,8,4)]:
        H=physical_graph(b,M,r)
        ev=np.linalg.eigvalsh(H)
        exact=exact_levels(b,M,r)
        assert len(exact)==len(ev)
        err=float(np.max(np.abs(ev-exact)))
        assert err<1e-10
        out.append(dict(b=b,M=M,r=r,L=b*M+r,dimension=len(ev),
                        complete_spectrum_error=err,
                        ground_energy=float(ev[0]),gap=float(ev[1]-ev[0])))
    return out


def entropy_from_blocks(blocks):
    eigen=[]
    maxdim=[0,0]
    for entries in blocks.values():
        left={a:i for i,a in enumerate(sorted({a for a,b,z in entries}))}
        right={b:i for i,b in enumerate(sorted({b for a,b,z in entries}))}
        mat=np.zeros((len(left),len(right)))
        for a,b,z in entries:mat[left[a],right[b]]=z
        eig=svdvals(mat,check_finite=False)**2
        eigen.extend(eig[eig>1e-15])
        maxdim=max(maxdim,list(mat.shape),key=lambda x:x[0]*x[1])
    eigen=np.asarray(eigen)
    assert abs(eigen.sum()-1)<1e-10
    return float(-np.sum(eigen*np.log(eigen))),maxdim


def marker_entropy(b,M):
    r=M//2;L=b*M+r;cut=L//2
    ks=2*np.pi*(np.arange(r)-(r-1)/2)/M
    blocks={};norm=0.;count=0
    # Cache each determinant for all physical first-marker positions.
    for longs in itertools.combinations(range(M),r):
        gaps=np.full(M,b,dtype=int);gaps[list(longs)]+=1
        relative=np.r_[0,np.cumsum(gaps)[:-1]]
        orbitals=np.exp(1j*np.outer(longs,ks))/np.sqrt(M)
        amplitude=np.sqrt(M/L)*abs(np.linalg.det(orbitals))
        for first in range(gaps[-1]):
            sites=relative+first
            leftbits=sum(1<<int(x) for x in sites if x<cut)
            rightbits=sum(1<<int(x-cut) for x in sites if x>=cut)
            number=leftbits.bit_count()
            blocks.setdefault(number,[]).append((leftbits,rightbits,amplitude))
            norm+=amplitude*amplitude;count+=1
    S,maxdim=entropy_from_blocks(blocks)
    assert abs(norm-1)<1e-10
    return dict(b=b,M=M,r=r,L=L,cut=cut,physical_dimension=count,
                entropy=S,norm_error=abs(norm-1),largest_schmidt_block=maxdim)


def xx_entropy(M):
    r=M//2;cut=M//2
    ks=2*np.pi*(np.arange(r)-(r-1)/2)/M
    C=np.exp(1j*np.outer(np.arange(cut),ks))/np.sqrt(M)
    eig=np.linalg.eigvalsh(C@C.conj().T)
    eig=np.clip(eig,1e-15,1-1e-15)
    return float(np.sum(-eig*np.log(eig)-(1-eig)*np.log(1-eig)))


def continuum_checks():
    # nu=1/2, with kappa=2/pi. Exactly compatible rings have M even.
    out=[]
    for b in [1,4,9]:
        Ms=np.array([40,48,64,80,96,128,160],dtype=float)
        L=(b+.5)*Ms;v=2*(b+.5);K=(b+.5)**-2
        E=2*Ms/np.pi-2/np.sin(np.pi/Ms)
        coef=np.linalg.lstsq(np.stack([np.ones_like(L),L**-2,L**-4],axis=1),E/L,rcond=None)[0]
        c=-6*coef[1]/(np.pi*v)
        # Current state: ground Slater set, physical P=2pi M/L.
        # Its exact energy shift is E_kin[cos(P/M)-1].
        M=160.;length=(b+.5)*M;P=2*np.pi/(b+.5)
        gap=(-2/np.sin(np.pi/M))*(np.cos(P/M)-1)
        K_est=gap*length/(2*np.pi*v)
        # Charge stiffness from neighboring M sectors at fixed L.
        def e_at(m):
            r=length-b*m
            return 2*m/np.pi-2*np.sin(np.pi*r/m)/np.sin(np.pi/m)
        second=e_at(M+1)+e_at(M-1)-2*e_at(M)
        invK_est=second*length/(np.pi*v)
        out.append(dict(b=b,velocity=v,K=K,casimir_c=float(c),
                        current_K_estimate=float(K_est),
                        charge_inverse_K_estimate=float(invK_est),
                        exact_inverse_K=1/K))
        assert abs(c-1)<1e-5
        assert abs(K_est/K-1)<1e-3
        assert abs(invK_est*K-1)<.02
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--max-M',type=int,default=16)
    args=parser.parse_args()
    entropies=[]
    for M in range(8,args.max_M+1,4):
        row=marker_entropy(4,M);row['xx_control_entropy']=xx_entropy(M)
        entropies.append(row)
        print(json.dumps(row),flush=True)
    slopes=[]
    for x,y in zip(entropies[:-1],entropies[1:]):
        slopes.append(dict(M_pair=[x['M'],y['M']],
                           marker_c=3*(y['entropy']-x['entropy'])/np.log(y['L']/x['L']),
                           xx_control_c=3*(y['xx_control_entropy']-x['xx_control_entropy'])/np.log(y['M']/x['M'])))
    out=dict(full_spectra=spectral_checks(),continuum=continuum_checks(),
             physical_entropies=entropies,half_chain_entropy_slopes=slopes,
             limitation='Projected marker model only, no microscopic finite-u CFT test.')
    with open('length_fluid_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()

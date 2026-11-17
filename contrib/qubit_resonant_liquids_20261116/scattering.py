#!/usr/bin/env python3
"""Analytic half-line scattering and sparse few-particle detuning checks."""
import json
from pathlib import Path
import time
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse.linalg import eigsh
from threadpoolctl import threadpool_limits
import leading as Ld

HERE=Path(__file__).resolve().parent


def halfline(q,delta,Q,n=2000):
    p=Ld.parameters(q);m,t=p['m'],p['t'];J=2*t*np.cos(Q/2)
    V=np.zeros(n)
    if delta==0:
        b=m;V[0]=-2*t;length=None
    elif q==1:
        b=4+delta;V[0]={1:10/3,2:14/3,3:0}[delta]
        length=b-m*J/(J+V[0])
    elif delta==1:
        b=7;length=b-m
    else:
        b=5;V[1]=8;v=V[1]/J;length=b-m*(1-2*v)/(1+2*v)
    ev=eigh_tridiagonal(V,-J*np.ones(n-1),select='i',select_range=(0,1),eigvals_only=True)
    predicted=-2*t-J*J/(2*t) if delta==0 and Q else None
    if predicted is not None: assert abs(ev[0]-predicted)<1e-8*t
    if delta: assert ev[0]>=-2*J-1e-9
    return dict(q=q,residue_difference=delta,coarse_total_momentum=Q,min_separation=b,
        continuum_bottom=-2*J,lowest_finite_halfline=float(ev[0]),bound_energy=predicted,
        edge_scattering_length=length,halfline_sites=n)


def detuned(q,L,N,gamma,phase):
    t=Ld.parameters(q)['t'];tau=t*gamma*np.exp(1j*phase)
    states=Ld.basis(L,q,N);_,H=Ld.matrix(L,q,states,shift=False,tau=tau)
    vals,vec=eigsh(H/t,k=min(3,len(states)-2),which='SA',tol=2e-11,ncv=min(60,len(states)))
    order=np.argsort(vals);E=float(vals[order[0]]);v=vec[:,order[0]]
    ks=2*np.pi*np.arange(L)/L
    band=-2*np.cos(Ld.parameters(q)['m']*ks)+2*np.real(tau/t*np.exp(-1j*ks))
    e1=float(min(band));res=float(np.linalg.norm(H/t@v-E*v))
    # Contact and cross-residue weights are diagnostics, not a fit to a CFT.
    pure=np.array([len({j%Ld.parameters(q)['m'] for j in Ld.bits(int(s),L)})==1 for s in states])
    return dict(q=q,L=L,N=N,gamma=gamma,phase=phase,dimension=len(states),E_over_t=E,E1_over_t=e1,
        binding_relative_to_N_free=E-N*e1,residual=res,pure_residue_weight=float(sum(abs(v[pure])**2)))


def main():
    start=time.time();half=[];det=[]
    with threadpool_limits(limits=1):
        for q in (1,3):
            for delta in range(Ld.parameters(q)['m']):
                for Q in (0.,.3,1.,2.):half.append(halfline(q,delta,Q))
        # Exploratory detunings: q3 real positive tau protects two valleys under
        # inversion; q1 phase pi/4 gives two equal minima only at leading order.
        for q,phase in ((1,np.pi/4),(3,0.)):
            for gamma in (.02,.1,.3):
                for L in (48,72,96):
                    r=detuned(q,L,2,gamma,phase);det.append(r);print(r,flush=True)
            for L in (24,36):
                r=detuned(q,L,3,.1,phase);det.append(r);print(r,flush=True)
    out=dict(halfline=half,nearest_hop_detunings=det,runtime_seconds=time.time()-start,
        interpretation='Few-particle binding only; no thermodynamic many-body phase or microscopic CFT is established.')
    (HERE/'scattering.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()

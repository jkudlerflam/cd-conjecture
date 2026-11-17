#!/usr/bin/env python3
"""Fixed microscopic test point and sixth-order detuning consistency checks.

This proposes a test point, not a certified two-component phase. All
coefficients are independent of L; no production ground-state run is made.
"""
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import leading as Ld
from higher_jets import series, A
from relative_detuning import minima


def wrap(k):
    return float((k+np.pi)%(2*np.pi)-np.pi)


def main():
    start=time.time(); checks=[]
    with threadpool_limits(limits=1):
        for q,phase in ((1,np.pi/4),(3,0.)):
            p=Ld.parameters(q); tau=.1*p['t']*np.exp(1j*phase)
            for L in (10,12,13):
                states,K,_=series(L,q,order=6,tau=tau,mu6=-.37*p['t'])
                _,H=Ld.matrix(L,q,states,shift=False,tau=tau)
                H=H+sparse.diags([.37*p['t']*int(s).bit_count() for s in states])
                error=A.maxabs(K[2]-H)
                assert error<1e-8
                checks.append(dict(q=q,L=L,sixth_error=error))
    # q1 is preferred here because q3's real hopping has a momentum-zero
    # opposite-valley pair channel, allowed by its exact microscopic symmetry.
    q=1; eps=1/1000; gamma=.1; phase=np.pi/4; nu=.05
    bands=minima(lambda k:-2*np.cos(4*k)+2*gamma*np.cos(k-phase))
    e0=bands[0][0]; ks=sorted(k for e,k in bands if abs(e-e0)<1e-9)
    assert len(ks)==2
    mu6=e0+nu; tau=gamma*np.exp(1j*phase)
    u=2*eps**2-10*eps**4+(50-mu6)*eps**6
    hop=-eps**2+4*eps**4+(-24+tau)*eps**6
    pair=eps**2
    h=A.bond(v=1,u=u,t=hop,p=pair,a=eps,b=-eps)
    assert np.max(abs(h-h.conj().T))==0
    channels=[]
    for a in range(-4,5):
        for b in range(-4,5):
            if 1<=abs(a)+abs(b)<=4 and a+b>0:
                channels.append(dict(net_valley_creation=[a,b],degree=abs(a)+abs(b),
                                     momentum=wrap(a*ks[0]+b*ks[1])))
    channels.sort(key=lambda x:abs(x['momentum']))
    out=dict(status='Unvalidated fixed microscopic test point; no CFT claim',
        q=q,epsilon=eps,gamma=gamma,phase=phase,nu=nu,mu6=mu6,
        leading_valleys=ks,leading_band_bottom=e0,leading_kinetic_scale=eps**6,
        microscopic=dict(v=1,a=eps,b=-eps,u=u,t=[hop.real,hop.imag],p=pair,f=0),
        bond_real=h.real.tolist(),bond_imag=h.imag.tolist(),
        low_degree_number_changing_momenta=channels,
        sixth_detuning_checks=checks,runtime_seconds=time.time()-start)
    Path(__file__).with_name('candidate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()

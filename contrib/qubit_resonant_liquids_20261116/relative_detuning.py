#!/usr/bin/env python3
"""Large half-line two-body check at the same/opposite valley thresholds.

Direct relative-coordinate reduction of leading.py, energies divided by t.
The cutoff is a relative-distance box, not the microscopic system size.
"""
import json
from pathlib import Path
import time
import numpy as np
from scipy.linalg import eig_banded
from scipy.optimize import minimize_scalar
from threadpoolctl import threadpool_limits
from leading import parameters


def minima(fn):
    xs=np.linspace(-np.pi,np.pi,4097)[:-1];ys=fn(xs);ds=2*np.pi/4096
    ii=np.flatnonzero((ys<=np.roll(ys,1))&(ys<=np.roll(ys,-1)))
    rs=[minimize_scalar(fn,bounds=(xs[i]-ds,xs[i]+ds),method='bounded',options={'xatol':1e-14}) for i in ii]
    return sorted([(float(r.fun),float(r.x)) for r in rs])


def solve(q,gamma,phase,K,rmax):
    p=parameters(q);m,t=p['m'],p['t']
    rr=[r for r in range(m,rmax+1) if not (q==3 and r==4)];pos={r:i for i,r in enumerate(rr)};n=len(rr)
    ab=np.zeros((m+1,n))
    pots={4:-2,5:10/3,6:14/3} if q==1 else {3:-2*t,8:8}
    for r,v in pots.items():ab[0,pos[r]]=v/t
    coeff={m:-2*np.cos(m*K/2),1:2*gamma*np.cos(K/2-phase)}
    for j,r in enumerate(rr):
        for dr,v in coeff.items():
            if r+dr in pos:ab[pos[r+dr]-j,j]=v
    ev,vec=eig_banded(ab,lower=True,select='i',select_range=(0,1))
    free=lambda k:sum(2*v*np.cos(d*k) for d,v in coeff.items())
    bottom=minima(free)[0][0]
    prob=abs(vec[:,0])**2
    return dict(q=q,gamma=gamma,phase=phase,K=K,rmax=rmax,
        continuum_bottom=bottom,E0=float(ev[0]),binding=float(ev[0]-bottom),
        short_distance_weight=float(sum(prob[np.array(rr)<=32])),mean_distance=float(np.dot(rr,prob)))


def main():
    start=time.time();rows=[]
    with threadpool_limits(limits=1):
        for q,phi in ((1,np.pi/4),(3,0.)):
            m=parameters(q)['m']
            for gamma in (.02,.1,.3):
                mins=minima(lambda k:-2*np.cos(m*k)+2*gamma*np.cos(k-phi))
                ks=[x for e,x in mins if abs(e-mins[0][0])<1e-8]
                assert len(ks)==2,(q,gamma,mins)
                for channel,K in (('same_valley',2*ks[0]),('opposite_valleys',ks[0]+ks[1])):
                    for rmax in (512,1024,2048):
                        row=solve(q,gamma,phi,K,rmax);row['channel']=channel
                        row['one_particle_minima']=ks;row['one_particle_energy']=mins[0][0]
                        rows.append(row);print(row,flush=True)
    Path(__file__).with_name('relative_detuning.json').write_text(json.dumps(dict(rows=rows,runtime_seconds=time.time()-start),indent=2)+'\n')


if __name__=='__main__':main()

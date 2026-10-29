"""Small periodic ED in the strict uniform-Z2 stratum. Run only in C67."""
import ed_core as ed
import numpy as np
import json
import argparse
from pathlib import Path
from scipy.optimize import brentq


def couplings(v):
    x,y,z,d,b=v
    J=np.diag([x,y,z]); J[0,1]=d; J[1,0]=-d
    return J,np.array([0.,0.,b])


def run(name,v,L,k=4):
    J,h=couplings(v)
    levels=[]; best=None; residual=0.; maxdim=0
    for p in [0,1]:
        data=ed.translation_orbits(L,p)
        for m in range(L):
            sec=ed.Sector(L,m,p,data)
            e,u,r,dim=ed.solve(sec,J,h,n=k,vector=True)
            levels.append(dict(p=p,m=m,e=e.tolist()))
            residual=max(residual,r); maxdim=max(maxdim,dim)
            if best is None or e[0]<best[0]-1e-9:
                best=(float(e[0]),p,m,sec.expand(u))
    e,p,m,psi=best
    S=ed.entropies(psi,L)
    cuts=np.arange(max(3,L//4),L//2+1)
    x=np.log(L/np.pi*np.sin(np.pi*cuts/L))
    c=float(3*np.polyfit(x,[S[str(l)] for l in cuts],1)[0])
    return dict(name=name,v=v,L=L,E0=e,p=p,m=m,S=S,c_chord=c,
                levels=levels,residual=residual,maxdim=maxdim)


def mass_tune(L,D,z=.25):
    # Ising x_epsilon/x_sigma = 8, restricted to m=0. This is a diagnostic,
    # not a thermodynamic definition once the chiral branch crosses below it.
    secs=[ed.Sector(L,0,p) for p in [0,1]]
    def f(g):
        J,h=couplings([-1.,-.5,z,D,-g])
        ee=ed.solve(secs[0],J,h,n=2)[0]
        eo=ed.solve(secs[1],J,h,n=1)[0][0]
        return float(ee[1]-ee[0]-8*(eo-ee[0]))
    g=brentq(f,1.90,2.08,xtol=2e-9)
    return g


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('mode',choices=['grid','junction','fixed','controls'])
    ap.add_argument('--sizes',nargs='+',type=int,default=[10,12,14,16])
    ap.add_argument('--names',nargs='*'); ap.add_argument('--Ds',nargs='+',type=float,default=[.1,.16,.20,.24])
    ap.add_argument('--out',required=True); a=ap.parse_args()
    points={}
    if a.mode=='grid':
        for name,r in json.loads(Path('grid_candidates.json').read_text()).items():
            points[name]=[r['v'][i] for i in [0,1,2,5,8]]
    elif a.mode=='controls':
        points={'ising':[-1.,0.,0.,0.,-1.], 'xx':[-1.,-1.,0.,0.,0.],
                'free_junction':[-1.,-.5,0.,.25,-1.5]}
    elif a.mode=='fixed':
        points={f'fixed_D{D:g}':[-1.,-.5,.25,D,-1.9891] for D in a.Ds}
    if a.names: points={k:v for k,v in points.items() if k in a.names}
    for L in a.sizes:
        jobs=points.copy()
        if a.mode=='junction':
            jobs={f'junction_D{D:g}':[-1.,-.5,.25,D,-mass_tune(L,D)] for D in a.Ds}
        for name,v in jobs.items():
            r=run(name,v,L)
            with open(a.out,'a') as f: f.write(json.dumps(r)+'\n')
            print(json.dumps({k:r[k] for k in ['name','L','v','E0','p','m','c_chord','residual']}),flush=True)


if __name__=='__main__': main()

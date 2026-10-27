"""Second-stage finite-momentum gap screen using previously computed grounds."""
import json, time
from d2 import density, sector_operator, solve, dump
import numpy as np
data=json.load(open('results/shortlist.json'))+json.load(open('results/global_shortlist.json'))
data += [r for r in json.load(open('results/boundary_scan.json')) if r.get('c_eff',0)>1.05]
out=[]
for p in data:
    t=time.time();h=density(np.array(p['J']),np.array(p['field'])); vs=[]
    for row in p['rows']:
        L=row['L']; v=[]
        for sign in (-1,1):
            H=sector_operator(L,h,(row['momentum']+sign)%L); e,_=solve(H,k=1)
            v.append(float((e[0]-row['E0'])*L/(2*np.pi)))
        vs.append(v)
    v1,v2=np.mean(vs,axis=1); L1,L2=[r['L'] for r in p['rows']]
    zq=1-np.log(v2/v1)/np.log(L2/L1) if min(v1,v2)>1e-8 else float('nan')
    base=min(r['c'] for r in p['rows'])-.6*p['c_drift']-.1*min(abs(p['z_gap']-1),5)
    score2=base-.4*min(abs(zq-1),5) if np.isfinite(zq) else -100
    r={k:p[k] for k in ('id','family','J','field','c_eff','c_drift','z_gap')}
    r.update(q1_velocities=vs,z_q1=float(zq),score2=float(score2),sizes=[L1,L2],seconds=time.time()-t)
    out.append(r);dump('results/velocity_screen.json',out)
    print(p['id'],'c',round(p['c_eff'],5),'zq',round(zq,4),'score',round(score2,4),flush=True)

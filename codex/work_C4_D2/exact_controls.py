"""Analytic ferromagnetic control and independent invariance checks."""
import json, math
from d2 import density, full_operator, quick, entropy_fit
import numpy as np
from scipy.special import gammaln

def logchoose(n,k):
    if k<0 or k>n: return -np.inf
    return gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
rows=[]
for L in range(8,23,2):
    S=[]
    for ell in range(1,L//2+1):
        logp=np.array([logchoose(ell,a)+logchoose(L-ell,L//2-a)-logchoose(L,L//2) for a in range(ell+1)])
        p=np.exp(logp); good=p>0; S.append(float(-p[good]@logp[good]))
    rows.append(dict(L=L,entropy=S,c_eff=entropy_fit(S,L)['c'],c_bulk=entropy_fit(S,L,max(2,L//4))['c'],
                     E0=-L,one_magnon_gap=4*(1-np.cos(2*np.pi/L)),c_MV=0.))
rng=np.random.default_rng(83); J=rng.normal(size=(3,3)); f=rng.normal(size=3)
Q,_=np.linalg.qr(rng.normal(size=(3,3))); Q[:,0]*=np.linalg.det(Q)
a=quick(J,f,6); b=quick(Q.T@J@Q,Q.T@f,6)
rotation=dict(energy_error=abs(a['E0']-b['E0']),entropy_error=float(np.max(abs(np.array(a['entropy'])-b['entropy']))))
H1=full_operator(6,density(J,f,1)); H2=full_operator(6,density(J,f,.5))
rotation['density_gauge_H_error']=float(np.linalg.norm((H1-H2).toarray()))
assert max(rotation.values())<1e-10
out=dict(ferromagnet=rows,invariance=rotation)
json.dump(out,open('results/exact_controls.json','w'),indent=2)
print(rotation)
for r in rows: print(r['L'],r['c_eff'],r['c_bulk'],r['one_magnon_gap'])

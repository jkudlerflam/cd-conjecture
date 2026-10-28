"""Independent exact-diagonalization check of the uniform contact commutator.
The analytic variational bound is in README.md. Numpy/scipy only.
"""
import itertools, json, math, time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

started=time.time()
rows=[]
for L in (6,8,10,12,14,16):
    N=L//2
    basis=[sum(1<<j for j in c) for c in itertools.combinations(range(L),N)]
    ix={x:i for i,x in enumerate(basis)}
    contact=np.array([sum(((x>>j)&1)*((x>>((j+1)%L))&1) for j in range(L)) for x in basis],float)
    rr=[];cc=[];vv=[]
    for a,x in enumerate(basis):
        for j in range(L):
            k=(j+1)%L
            if ((x>>j)&1)!=((x>>k)&1):
                rr.append(ix[x^(1<<j)^(1<<k)]);cc.append(a);vv.append(-1.)
    T=coo_matrix((vv,(rr,cc)),shape=(len(basis),)*2).tocsr()
    en,vec=eigsh(T,k=1,which='SA',tol=1e-12)
    e0=float(en[0]);psi=vec[:,0]
    Cp=contact*psi
    B=2*np.vdot(Cp,T@Cp-e0*Cp).real
    exact=0.
    # Filled sea is antiperiodic for N even, periodic for N odd.
    ks=2*np.pi/L*(np.arange(N)-(N-1)/2)
    def corr(r):return np.exp(1j*ks*r).sum()/L
    mat=np.array([[corr(a-b) for b in (-1,2,1)] for a in (-1,2,0)])
    pair1=.5*corr(-1)-corr(-2)*corr(1)
    pair2=.5*corr(-1)-corr(1)*corr(-2)
    wick=2*(pair1+pair2-2*np.linalg.det(mat)).real
    ev=[]
    for v in (-.04,-.02,.02,.04):
        eig=eigsh(T+diags(v*contact),k=1,which='SA',tol=1e-12,v0=psi)[0][0]
        delta=float(eig-e0-v*np.dot(psi*psi,contact))
        guaranteed=-v*v*(B/L)**2*L/(360*(1+abs(v)))
        ev.append(dict(v=v,relaxation=delta,relaxation_per_physical_site_over_v2=delta/L/v**2,variational_upper_bound=guaranteed,passes=bool(delta<=guaranteed+1e-10)))
    rows.append(dict(L=L,N=N,dimension=len(basis),energy=e0,contact_per_material_site=float(np.dot(psi*psi,contact)/N),B_per_physical_site=float(B/L),wick_B_per_physical_site=float(wick),wick_error=float(abs(wick-B/L)),residual=float(np.linalg.norm(T@psi-e0*psi)),checks=ev))
result=dict(description='Uniform scalar contact ground-energy curvature and local-commutator bound; not a c=24 simulation.',thermodynamic_B_per_physical_site=1/math.pi+16/(9*math.pi**3),seconds=time.time()-started,rows=rows)
p=Path(__file__).with_name('uniform_contact_results.json');p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

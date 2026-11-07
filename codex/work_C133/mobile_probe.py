"""Full d=6 mobile toy at short lengths, plus exact legal-sector comparisons.
The imported source is read only and bytecode writing is disabled by -B.
This is a finite-size control, not an asymptotic A measurement.
"""
import importlib.util,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scan import ground
spec=importlib.util.spec_from_file_location('c133_mobile_source',Path('../../contrib/mobile_cells/ed.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def normalized_width(h,o):
    x=h+np.kron(o,np.eye(6))-np.kron(np.eye(6),o)
    e=np.linalg.eigvalsh(x);return (e[-1]-e[0])/2

def bracket_M(h):
    # A full gauge lower bound from equal one-site marginals of a maximally
    # entangled vector and a product |aa>. Their expectation gap is gauge invariant.
    d=6;phi=np.eye(d).ravel()/np.sqrt(d)
    trial=[np.eye(d*d)[:,d*a+a] for a in range(d)]+[phi]
    energies=[float(v@h@v) for v in trial]
    lower=(max(energies)-min(energies))/2
    # Diagonal-gauge minimization is only an upper bound, and is labeled as such.
    fun=lambda z:normalized_width(h,np.diag(np.r_[z,0.]))
    opt=minimize(fun,np.zeros(5),method='Powell',options={'xtol':1e-10,'ftol':1e-11,'maxiter':1000})
    return lower,float(opt.fun)

rows=[]
for u in [.1,.3,1.]:
    p=dict(t=1.,u=u,mu=-1.,J=.5,Jz=0.,K=10.)
    h=m.dense_bond(**p);lo,hi=bracket_M(h);energies={}
    for L in [3,4,6]:
        e,res,n=ground(h,6,L);energies[L]=e
        legal=[0.]
        for W in range(2,L+1,2):
            for nup in range(W+1):
                r=m.Ring(L,W,nup,t=1.,u=u,mu=-1.,J=.5,Jz=0.)
                H=r.sparse()
                if r.dim<=2:ee=np.linalg.eigvalsh(H.toarray())[0]
                else:
                    from scipy.sparse.linalg import eigsh
                    ee=eigsh(H,k=1,which='SA',return_eigenvectors=False,tol=1e-11)[0]
                legal.append(float(ee))
        assert abs(e-min(legal))<1e-8
        rows.append(dict(u=u,L=L,E=e,M_lower=lo,M_upper=hi,residual=res,legal_error=abs(e-min(legal))))
    Qnum=3*(energies[6]-2*energies[3]);print(u,lo,hi,'Q interval',sorted([Qnum/lo,Qnum/hi]),flush=True)
Path('mobile_results.json').write_text(json.dumps(rows,indent=2)+'\n')
text=['# Mobile toy short-ring check','','[numerical] Full d=6 bond from contrib/mobile_cells/ed.py, with t=1, mu=-1, J=0.5, Jz=0, K=10, and varying marker hop u. All legal wall-number and spin sectors were checked against unrestricted ED at L=3,4,6. They agree. Gauge optimization is bracketed, not asserted solved. The lower bound uses gauge-invariant expectations in equal-marginal states, the upper bound minimizes over diagonal coboundaries.','','[numerical, limitation] These short rings are not certified vacuum-compatible along the incommensurate charge fluid. Therefore Q_3 is recorded as a raw finite-ring diagnostic and no A estimate is assigned. Increasing marker mobility adds at most the one length mode in this toy. It has no parameter encoding an unbounded payload. Larger sector-resolved work is needed to infer its microscopic A/M.','','| u | M interval | raw Q_3 interval |','| --: | :-- | :-- |']
for u in [.1,.3,1.]:
    rs=[r for r in rows if r['u']==u];by={r['L']:r for r in rs};lo=rs[0]['M_lower'];hi=rs[0]['M_upper'];num=3*(by[6]['E']-2*by[3]['E']);qs=sorted([num/lo,num/hi]);text.append(f'| {u} | [{lo:.9g}, {hi:.9g}] | [{qs[0]:.9g}, {qs[1]:.9g}] |')
text+=['','[proved] For any vector with identical one-site marginals, the expectation of o tensor I minus I tensor o vanishes. Two such vectors therefore certify M >= half their expectation difference. This lower certificate applies to every Hermitian gauge, even though the numerical upper search only used diagonal gauges.']
Path('MOBILE.md').write_text('\n'.join(text)+'\n')

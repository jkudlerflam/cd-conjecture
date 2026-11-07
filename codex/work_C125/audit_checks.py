"""Independent checks of the algebraic claims, with explicit failure witnesses."""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction
from math import factorial
import json
import numpy as np
from scipy.special import logsumexp
from scipy.linalg import expm
from checks import ring_spectrum, transfer, bond
ROOT=Path(__file__).resolve().parent

def rising(a,r):
    out=1
    for i in range(r):out*=a+i
    return out

def racah(k,l,J):
    return sum((Fraction(rising(-l,r)*rising(l+1,r)*rising(-J,r)*rising(J+1,r),
                 rising(1,r)*rising(-k,r)*rising(k+2,r)*factorial(r))
                for r in range(min(l,J)+1)),Fraction(0))/(k+1)

recoupling=[]
for k in range(1,13):
    d=k+1;H=[sum((Fraction(1,i) for i in range(1,J+1)),Fraction(0)) for J in range(k+1)]
    for l in range(k+1):
        r=[racah(k,l,J) for J in range(k+1)]
        for J in range(k+1):
            a=Fraction((J+1)*(d*d-(J+1)**2),2)
            am=Fraction(J*(d*d-J*J),2)
            Dr=(a*(r[J]-r[J+1]) if J<k else 0)+(am*(r[J]-r[J-1]) if J else 0)
            assert Dr==l*(l+1)*(2*J+1)*r[J]
        lam=-sum((2*J+1)*r[J]*H[J] for J in range(k+1))
        expected=Fraction(d,2*l*(l+1)) if l else -d*H[k]+Fraction(k,2)
        assert lam==expected
    recoupling.append(dict(k=k,exact_recurrence=True,exact_coefficients=True))

# Momentum-resolved compact boson, j=2 Delta. Aggregate degeneracies, never allocate states.
jmax=24;part=[0]*(jmax//2+1);part[0]=1
for k in range(1,len(part)):
    for N in range(k,len(part)):part[N]+=part[N-k]
one=defaultdict(int)
for m in range(-5,6):
 for w in range(-5,6):
  for N in range(len(part)):
   for Nb in range(len(part)):
    j=m*m+w*w+2*(N+Nb)
    if j<=jmax:one[j,m*w+N-Nb]+=part[N]*part[Nb]
coef={(0,0):1};rows=[]
for n in range(1,5):
    new=defaultdict(int)
    for (j,p),c in coef.items():
        for (jj,pp),cc in one.items():
            if j+jj<=jmax:new[j+jj,p+pp]+=c*cc
    coef=new
    if n not in [1,2,4]:continue
    shells=defaultdict(int)
    for (j,p),c in coef.items():shells[j]+=c
    for beta in [2.,4.,8.]:
        rank=0;last=0
        for j in sorted(shells):
            if np.log(rank+shells[j])>beta/2:break
            rank+=shells[j];last=j
        energies=[];momenta=[];multip=[]
        for (j,p),c in coef.items():
            if j<=last:
                energies.append(np.pi*j/beta);momenta.append(2*np.pi*p/beta);multip.append(c)
        e=np.array(energies);p=np.array(momenta);lm=np.log(np.array(multip,float))
        for L in [1.,4.,16.]:
            g=logsumexp(lm-L*e);g2=logsumexp(lm-2*L*e)
            w=np.exp(lm-L*e-g);w2=np.exp(lm-2*L*e-g2)
            derivative=2*(np.dot(w2,e)-np.dot(w,e))
            for f in [0.,.5,.99]:
                a=f*L/(4*np.e);z=logsumexp(lm-L*e-a*p);wt=np.exp(lm-L*e-a*p-z)
                X=np.vstack([e,p]);mean=X@wt;Y=X-mean[:,None];cov=(Y*wt)@Y.T
                rows.append(dict(n=n,beta=beta,L=L,a=a,rank=rank,
                    cone_margin=float(min(e-abs(p)/(4*np.e))),
                    hessian_min_eigenvalue=float(min(np.linalg.eigvalsh(cov))),
                    R=float(2*g-g2),R_derivative=float(derivative)))
                assert min(np.linalg.eigvalsh(cov))>=-1e-10
                assert min(e-abs(p)/(4*np.e))>=-1e-12
                assert derivative<=1e-12
                assert -1e-12<=2*g-g2<=np.log(rank)+1e-12

# Test monotonicity away from the two-site ring convention for attractive XXZ.
xxz=[]
for delta in [-.9,-.5]:
    spec={L:ring_spectrum(2,L,delta) for L in [4,6,8,12]}
    best=(-np.inf,None)
    for beta in np.geomspace(.1,256,100):
        g={L:logsumexp(-beta*e) for L,e in spec.items()}
        R4=2*g[4]-g[8];R6=2*g[6]-g[12]
        if R6-R4>best[0]:best=(R6-R4,dict(delta=delta,beta=float(beta),R4=float(R4),R6=float(R6),increase=float(R6-R4)))
    xxz.append(best[1])

# Independent contraction of the original length-four Trotter trace.
trotter=[]
for name,d,delta in [('XX',2,0.),('XXZ_m09',2,-.9),('Heisenberg',2,1.),('SU3',3,0.)]:
    h=bond(delta) if d==2 else np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9)
    L=4;powers=d**np.arange(L);Hs=[np.zeros((d**L,d**L)) for _ in range(2)]
    for state in range(d**L):
        digits=(state//powers)%d
        for site in range(L):
            j=(site+1)%L;a,b=int(digits[site]),int(digits[j])
            for aa in range(d):
                for bb in range(d):
                    target=state+(aa-a)*powers[site]+(bb-b)*powers[j]
                    Hs[site%2][target,state]+=h[d*aa+bb,d*a+b]
    beta=2.;m=2
    z1=np.trace(np.linalg.matrix_power(expm(-beta*Hs[0]/m)@expm(-beta*Hs[1]/m),m))
    U=transfer(d,h,beta,m);z2=np.trace(U@U).real
    residual=abs(z1-z2)/z1
    assert residual<1e-11
    trotter.append(dict(model=name,beta=beta,m=m,relative_residual=float(residual)))

base=json.loads((ROOT/'checks.json').read_text())
for r in base['matsubara']:
    V=2*np.e;omega=r['omega'];alpha=omega/V if omega<=V else np.log(omega/2)
    assert r['rate']>=alpha-1e-12
assert base['xx_ED_max_error']<1e-10
assert max(abs(r['error']) for r in base['xx_crossed_products'])<1e-7
assert all(r['bound_margin']>=-1e-10 for r in base['xx_large'])
assert all(r['trivial_margin']>=-1e-10 and r['gluing_margin']>=-1e-10 for r in base['finite_trace_checks'])
assert all(r['max_cA_increase']<0 and r['min_theta_slope_increase']>0 for r in base['bethe_sym2'])
assert not any(r['enumeration_exhausted'] for r in base['countermodel_finite'])
assert base['ring_curvatures']['XXZ_m09']['min_theta_slope_increase']<0
assert base['doubling_ground']['SU3']['4']<0
result=dict(trotter_contraction=trotter,exact_recoupling=recoupling,countermodel_covariances=rows,attractive_R_monotonicity=xxz)
(ROOT/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(recurrences_checked=len(recoupling),countermodel_checks=len(rows),attractive_R_monotonicity=xxz),indent=2))

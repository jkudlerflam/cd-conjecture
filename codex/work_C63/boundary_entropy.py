"""Schmidt eigenvalues of the exact filtered-Dicke ground state at u=r.
No many-body diagonalization. All coefficients use numpy/scipy.
"""
import json
import numpy as np
from scipy.special import gammaln
from pathlib import Path


def coeff(ne,no,q):
    def binom(n,x):
        k=np.arange(n+1)
        return np.exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+k*np.log(x)+(n-k)*np.log1p(-x))
    return np.convolve(binom(ne,q),binom(no,1-q))


def S(L,ell,b,start=0):
    r=np.sqrt(1+b*b/4);w=(2*r+b)/(2*r-b);q=np.sqrt(w)/(1+np.sqrt(w))
    ne=(ell+1-start)//2;no=ell-ne
    ca=coeff(ne,no,q);cb=coeff(L//2-ne,L//2-no,q)
    m=np.arange(ell+1);idx=L//2-m;valid=(idx>=0)&(idx<len(cb))
    lam=ca[valid]*cb[idx[valid]];lam/=np.sum(lam);lam=lam[lam>1e-300]
    return float(-np.sum(lam*np.log(lam)))

out=[]
for b in [0.,2.,4.]:
    for L in [16,32,64,128,256,512]:
        ell=np.arange(L//4,L//2+1)
        ss=np.array([.5*(S(L,int(l),b)+S(L,int(l),b,1)) for l in ell])
        cch=3*np.polyfit(np.log(np.sin(np.pi*ell/L)),ss,1)[0]
        cs=3*(S(L,L//2,b)-S(L//2,L//4,b))/np.log(2)
        r=np.sqrt(1+b*b/4);gap=4*r-4*np.sqrt(r*r-np.sin(2*np.pi/L)**2)
        row=dict(b=b,L=L,c_chord=float(cch),c_size=float(cs),S_half=S(L,L//2,b),exact_defect_gap=float(gap),v_defect_est=float(L*gap/(2*np.pi)))
        out.append(row);print(json.dumps(row))
Path('boundary_entropy.json').write_text(json.dumps(out,indent=2)+'\n')

# Independent finite-ring verification of the filtered state and Schmidt formula.
from ed_staggered import sector,entropy_profile
checks=[]
for b in [0.,2.,4.]:
    L=8;r=np.sqrt(1+b*b/4);H,states=sector(L,L//2,r,b)
    even=np.sum((states[:,None]>>np.arange(0,L,2))&1,axis=1)
    ratio=np.sqrt((2*r+b)/(2*r-b));psi=(-ratio)**even;psi/=np.linalg.norm(psi)
    residual=float(np.linalg.norm(H@psi+r*L*psi))
    ent,_,_=entropy_profile(psi,states,L)
    error=max(abs(ent[l-1]-.5*(S(L,l,b)+S(L,l,b,1))) for l in range(1,L//2+1))
    assert max(residual,error)<1e-10
    checks.append(dict(b=b,L=L,ground_state_residual=residual,entropy_error=error))
Path('boundary_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print('Boundary state checks',checks)

"""Rational selector inequalities and small checks of the uniform bound."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from fractions import Fraction as F
import json
import numpy as np
from scipy.sparse.linalg import eigsh
import bond15,bond18

def atan_bounds(x,terms=60):
    s=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    t=(-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return min(s,s+t),max(s,s+t)
def cos_point(x,terms=30):
    s=F(1);t=F(1)
    for j in range(1,terms):t=-t*x*x/F((2*j-1)*(2*j));s+=t
    nxt=-t*x*x/F((2*terms-1)*(2*terms))
    return min(s,s+nxt),max(s,s+nxt)
class I:
    def __init__(self,l,u=None):self.l=F(l);self.u=F(l if u is None else u)
    def __add__(self,x):
        if not isinstance(x,I):x=I(x)
        return I(self.l+x.l,self.u+x.u)
    __radd__=__add__
    def __neg__(self):return I(-self.u,-self.l)
    def __sub__(self,x):return self+-x
    def __rsub__(self,x):return -self+x
    def __mul__(self,x):
        assert x>=0
        return I(self.l*x,self.u*x)
    __rmul__=__mul__
    def __truediv__(self,x):return self*F(1,x)
def main():
    a,b=atan_bounds(F(1,5));c,d=atan_bounds(F(1,239))
    pi=I(16*a-4*d,16*b-4*c)
    def cospi(den,num=1):return I(cos_point(pi.u*num/den)[0],cos_point(pi.l*num/den)[1])
    aa=lambda b:2*(1-cospi(b+1))
    dd=lambda b:(b+1)*aa(b)-b*aa(b+1)
    mu=2-(dd(16)+dd(15))/2
    ee=lambda b:(mu-2*cospi(b+1))/b
    sep=[ee(15)-ee(16),ee(17)-ee(16)]
    assert min(x.l for x in sep)>F(1,1<<16)
    assert ee(16).u<-F(1,1<<8)
    gap=2*(cospi(17)-cospi(17,2));assert gap.l>F(1,16)
    assert (2-2*cospi(17)).l>F(1,64)
    out=dict(selector_density_gap_lower=float(min(x.l for x in sep)),selected_density_upper=float(ee(16).u),
             mu_sel=[float(mu.l),float(mu.u)],vacuum_gap_lower=float(gap.l),
             rational_acceptance=True,small_checks=[])
    for n in [1,2]:
        p=bond18.params(n,K=100.,g=.03)
        z=bond15.cluster(p);H=z['HH'];bad=z['bad']
        # Keep precisely the matrix entries preserving every bad-pair
        # projector. For these rules, a vacuum bad move changes N_bad by 1,
        # and every other off-diagonal move preserves all bad projectors.
        coo=H.tocoo();keep=bad[coo.row]==bad[coo.col]
        import scipy.sparse as sp
        Hg=sp.coo_matrix((coo.data[keep],(coo.row[keep],coo.col[keep])),shape=H.shape).tocsc()
        ideal=2*(p['musel']-2*np.cos(np.pi/(n+2)))
        # H_good >= ideal-4 N_bad for the pair; K remains in Hg.
        test=Hg-sp.diags((p['K']-4)*bad)-ideal*sp.eye(H.shape[0])
        low=float(eigsh(test,k=1,sigma=-.1,which='LM',ncv=40,tol=1e-10,return_eigenvectors=False)[0])
        assert low>-1e-8
        actual=float(eigsh(H,k=1,sigma=ideal-.1,which='LM',return_eigenvectors=False)[0])
        lower=ideal-2*max(0,n-1)/(p['K']-4)
        assert actual>lower-1e-9
        out['small_checks'].append(dict(n=n,dimension=H.shape[0],good_form_min=low,energy=actual,uniform_lower_bound=lower))
    print(json.dumps(out),flush=True);json.dump(out,open('selector_certificate.json','w'),indent=2)
if __name__=='__main__':main()

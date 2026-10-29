"""C85 rational selector at b=15, reusing C84 alternating enclosures."""
from fractions import Fraction as F
import json
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
    mu=2-(dd(15)+dd(14))/2
    ee=lambda b:(mu-2*cospi(b+1))/b
    sep=[ee(14)-ee(15),ee(16)-ee(15)]
    delta=min(x.l for x in sep)
    gap=2*(cospi(16)-cospi(16,2))
    assert delta>F(1,1<<16) and ee(15).u<-F(1,1<<8)
    assert gap.l>F(1,16)
    assert delta-F(1,(1<<2560)-4)>F(1,1<<17)
    assert F(1,1<<17)>F(8,1<<608)
    out=dict(cell_length=15,mu_sel=[float(mu.l),float(mu.u)],
      selector_density_gap_lower=float(delta),selected_density_upper=float(ee(15).u),
      vacuum_gap_lower=float(gap.l),F3_verified=True,rational_acceptance=True)
    def dyadic(x,upper=False):
        y=x*(1<<128);k=y.numerator//y.denominator
        return F(k+int(upper and y.denominator!=1),1<<128)
    dl=dyadic(delta);ml=dyadic(mu.l);mu_hi=dyadic(mu.u,True)
    exact=dict(delta_lower=[str(dl.numerator),str(dl.denominator)],
      mu_sel=[[str(x.numerator),str(x.denominator)] for x in [ml,mu_hi]])
    json.dump(out,open('selector_certificate.json','w'),indent=2)
    json.dump(exact,open('selector_exact.json','w'),indent=2)
    print(json.dumps(out))
if __name__=='__main__':main()

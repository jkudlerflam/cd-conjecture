"""Rational interval certificate for the full C88 response obstruction."""
from fractions import Fraction as F
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    p=[x*y for x in a for y in b];return min(p),max(p)
def div(a,b):
    assert b[0]>0
    return mul(a,(1/b[1],1/b[0]))
def const(x):return (F(x),F(x))
def atan_inv(q,m=40):
    s=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(m)),F(0))
    e=F(1,(2*m+1)*q**(2*m+1))
    return s,s+e # m even, partial sum below the alternating limit
pi=sub(mul(const(16),atan_inv(5)),mul(const(4),atan_inv(239)))
def sine(x,m=36):
    term=x;total=x
    for k in range(1,m):
        term=div(neg(mul(term,mul(x,x))),const((2*k)*(2*k+1)))
        total=add(total,term)
    # Taylor remainder bounded by |x|^(2m+1)/(2m+1)!.
    import math
    rem=max(abs(x[0]),abs(x[1]))**(2*m+1)/math.factorial(2*m+1)
    return total[0]-rem,total[1]+rem

def energy(n,v):
    if v==0:return const(0)
    return mul(const(-2),div(sine(div(mul(const(n),pi),const(n+v))),sine(div(pi,const(n+v)))))
def contact(n,v):
    l=n+v
    sn=sine(div(mul(const(n),pi),const(l)));sd=sine(div(pi,const(l)))
    return sub(const(F(n,l)),div(mul(sn,sn),mul(const(n*l),mul(sd,sd))))
def out(a):return [float(x) for x in a]

def main():
    n=15;v=15
    ev=energy(n,v);ew=energy(n,v+1);s=sub(ev,ew)
    baseline=add(ev,mul(const(v),s))
    margins=[]
    for w in range(40):
        if w in (v,v+1):continue
        delta=sub(add(energy(n,w),mul(const(w),s)),baseline)
        assert delta[0]>0
        margins.append((w,delta))
    tail=sub(add(const(-2*n),mul(const(40),s)),baseline)
    assert tail[0]>0 and s[0]>F(1,2) and s[1]<F(3,4)
    da=sub(contact(n,v),contact(n,v+1));assert da[0]>0
    # Full legal trace: five length-3 cells have mean D energy at least
    # 5*(6/13)*1000. The single length-15 cell has norm below 1016.
    moment_lower=F(30000,13)**2-F(1016)**2
    slope=mul(div(mul(da,da),const(4)),const(moment_lower))
    assert slope[0]>0
    result=dict(tag='proved_by_rational_interval_arithmetic',N=n,V_pair=[v,v+1],
        tilt_interval=out(s),contact_difference_interval=out(da),
        other_vacancy_margin_lower=float(min(d[0] for _,d in margins)),
        tail_from_V=40,tail_margin_lower=float(tail[0]),
        moment_difference_strict_lower=float(moment_lower),
        response_T_times_difference_strict_lower=float(slope[0]),
        selected_lengths=[15],wrong_lengths=[3]*5,W=15)
    (HERE/'crossing_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()

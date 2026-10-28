"""Explicit all-n rational schedule, evaluated only through logarithms."""
import math,json
from fractions import Fraction
from structured import visibility

def lg(x):
    if isinstance(x,Fraction):return math.log10(x.numerator)-math.log10(x.denominator)
    return math.log10(x)

def schedule(n):
    det=Fraction(1)
    for i in range(n):
        for j in range(i+1,n):det*=4*(j-i)**2
        for j in range(n):det/=2*i+2*j+3
    trace=sum((Fraction(1,4*i+3) for i in range(n)),Fraction(0))
    q=det/trace**(n-1);B=32*n*4**n
    ell=max(2,math.ceil((math.log2(4*n*B*B)-(lg(q)/math.log10(2)))/4))
    while Fraction(4*n*B*B,1<<(4*ell))>q:ell+=1
    L=q/Fraction(32*math.factorial(2*n-1)**2*(1<<(ell*(4*n-1))))
    Dfrac=4*n*(n+2)/L
    D=(Dfrac.numerator+Dfrac.denominator-1)//Dfrac.denominator
    s=L/2
    lb=-(n+1)*lg(D+6)
    logg=lb+lg(s)-lg(16*(n+1))
    loggamma=math.log10(visibility(n))+2*logg+2*lb+lg(s)
    return dict(n=n,log10_q=lg(q),t0_power_of_two=-ell,log10_gram_lower=lg(L),
                log10_D=lg(D),log10_beta_lower=lb,log10_g=logg,
                log10_herm_margin_lower=loggamma,log10_field=loggamma-math.log10(2))
if __name__=='__main__':
    out=[schedule(n) for n in [1,2,3,4,8,12,16,23,32,48]]
    json.dump(out,open('schedule_results.json','w'),indent=2)
    for r in out:print(json.dumps(r))

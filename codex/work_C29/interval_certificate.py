"""Rational interval certificate for the four reference Fermi points."""
import json
from fractions import Fraction as Q
from itertools import combinations
from checks import rational_weights


def atan_bounds(x,N=32):
    terms=[(-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(N)]
    s=sum(terms);nextterm=(-1)**N*x**(2*N+1)/Q(2*N+1)
    return min(s,s+nextterm),max(s,s+nextterm)


def trig_bounds(x,sine=False,N=16):
    term=x if sine else Q(1);s=term
    for k in range(1,N):
        term*=-x*x/Q((2*k)*(2*k+1) if sine else (2*k-1)*(2*k));s+=term
    nxt=term*(-x*x/Q((2*N)*(2*N+1) if sine else (2*N-1)*(2*N)))
    return min(s,s+nxt),max(s,s+nxt)


def tan_bounds(x):
    s=trig_bounds(x,True);c=trig_bounds(x)
    assert c[0]>0
    return s[0]/c[1],s[1]/c[0]


def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def times(k,a):return (k*a[0],k*a[1]) if k>=0 else (k*a[1],k*a[0])


def main():
    aa=atan_bounds(Q(1,5));bb=atan_bounds(Q(1,239));pi=sub(times(16,aa),times(4,bb))
    f0,f1=rational_weights(3);h0,h1=rational_weights(5)
    xs=[(Q('0.8241'),Q('0.8242')),(Q('0.8396'),Q('0.8397'))]
    for x,y in zip(xs,(64*h0/f0,64*h1/f1)):
        assert tan_bounds(x[0])[1]<y<tan_bounds(x[1])[0]
    ks=[(Q('1.8234'),Q('1.8235')),(Q('1.2132'),Q('1.2133'))]
    for k,y in zip(ks,(Q(-1,4),Q(7,20))):
        assert trig_bounds(k[0])[0]>y>trig_bounds(k[1])[1]
    points=[sub(add(pi,xs[0]),ks[0]),sub(add(xs[0],ks[0]),pi),sub(xs[1],ks[1]),add(xs[1],ks[1])]
    def separation(a,b):
        diff=sub(a,b)
        # Integer choice alone uses a float; the bound and its sign are rational.
        winding=round(float((diff[0]+diff[1])/(2*(pi[0]+pi[1]))))
        reduced=sub(diff,times(2*winding,pi))
        assert reduced[0]>-pi[0] and reduced[1]<pi[0]
        assert reduced[0]>0 or reduced[1]<0
        return min(abs(reduced[0]),abs(reduced[1]))
    bil=min(separation(points[i],points[j]) for i,j in combinations(range(4),2))
    pairs=[add(points[i],points[j]) for i,j in combinations(range(4),2)]
    quart=min(separation(pairs[i],pairs[j]) for i,j in combinations(range(6),2))
    print('rational lower bounds',float(bil),float(quart))
    assert bil>Q('0.088') and quart>Q('0.030')
    out=dict(status='proved by rational alternating-series bounds',point_intervals=[[float(x),float(y)] for x,y in points],bilinear_lower_bound='0.088',distinct_pair_sum_lower_bound='0.030')
    print(json.dumps(out,indent=2))
    with open('interval_certificate_results.json','w') as f:json.dump(out,f,indent=2)

if __name__=='__main__':main()

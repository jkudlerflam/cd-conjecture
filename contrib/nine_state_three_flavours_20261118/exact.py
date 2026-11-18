"""Rational carrier inverse, cubic root certificate, and a literal control box."""
from fractions import Fraction as F
from pathlib import Path
import itertools,json
import numpy as np
from physical import CONTROL_NUMERATORS,CONTROL_DENOMINATOR
from exact_symbol import inverse


class Interval:
    def __init__(self,lo,hi=None):self.lo=F(lo);self.hi=F(lo if hi is None else hi)
    def __add__(self,b):
        b=b if isinstance(b,Interval) else Interval(b)
        return Interval(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return Interval(-self.hi,-self.lo)
    def __sub__(self,b):return self+-b
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):
        b=b if isinstance(b,Interval) else Interval(b)
        z=[a*c for a in (self.lo,self.hi) for c in (b.lo,b.hi)]
        return Interval(min(z),max(z))
    __rmul__=__mul__


def exact_T():
    n=3;s=5;H=2*np.eye(s,dtype=int)-np.eye(s,k=1,dtype=int)-np.eye(s,k=-1,dtype=int)
    Q=np.kron(H,np.eye(s,dtype=int))+np.kron(np.eye(s,dtype=int),H)
    R=inverse(Q);assert np.array_equal(Q@R,np.eye(s*s,dtype=int))
    x=R[:,(s-1)*s].reshape(s,s);U=np.full((s*s,n),F(0),dtype=object);V=U.copy()
    for i in range(n):U[i*s:(i+1)*s,i]=x[i,:];V[i::s,i]=x[:,i]
    return -2*(U.T@R@V)


def exact_M(controls):
    mu,eta,b,L,R=controls
    return np.array([[mu+F(11,12)*(eta+b),F(2,3)*L+F(1,12)*R,F(0)],
                     [F(2,3)*L+F(1,12)*R,mu+F(2,3)*eta,F(1,3)*(L+R)],
                     [F(0),F(1,3)*(L+R),mu+F(1,3)*eta]],dtype=object)


def det3(X):
    result=0
    for p in itertools.permutations(range(3)):
        sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        result+=sign*X[0,p[0]]*X[1,p[1]]*X[2,p[2]]
    return result


def determinant_at(T,M,x):
    S=(T+T.T)/2;A=(T-T.T)/2;a=[A[1,2],-A[0,2],A[0,1]]
    X=M+2*x*S
    return det3(X)-4*(1-x*x)*sum(a[i]*X[i,j]*a[j] for i in range(3) for j in range(3))


def certificate():
    T=exact_T();controls=[F(int(x),CONTROL_DENOMINATOR) for x in CONTROL_NUMERATORS];M=exact_M(controls)
    y={i:determinant_at(T,M,F(i)) for i in (0,1,-1,2)}
    p0=y[0];p2=(y[1]+y[-1]-2*p0)/2;ac=(y[1]-y[-1])/2
    p3=(y[2]-p0-4*p2-2*ac)/6;p1=ac-p3;poly=[p0,p1,p2,p3]
    assert p3==8*det3(T) and p3!=0
    # Three disjoint sign-changing intervals exhaust a nonzero cubic.
    brackets=[(F(-7,10),F(-13,20)),(F(9,50),F(23,100)),(F(9,20),F(1,2))]
    endpoint_values=[[determinant_at(T,M,x) for x in interval] for interval in brackets]
    assert all(a*b<0 for a,b in endpoint_values)
    derivative_bounds=[]
    for lo,hi in brackets:
        x=Interval(lo,hi);d=p1+2*p2*x+3*p3*x*x
        assert d.lo*d.hi>0
        derivative_bounds.append(min(abs(d.lo),abs(d.hi)))
    norm_bound=max(sum(abs(x) for x in row) for row in M)+max(sum(abs(x) for x in row) for row in T)+max(sum(abs(x) for x in row) for row in T.T)
    velocity_lower=min(derivative_bounds)/(2*norm_bound**2)
    # The same six bracket signs for every independent physical control in a box.
    box_width=None
    for exponent in range(6,15):
        width=F(1,10**exponent);MI=exact_M([Interval(c-width,c+width) for c in controls])
        signs=[]
        for pair,values in zip(brackets,endpoint_values):
            for x,value in zip(pair,values):
                bound=determinant_at(T,MI,x)
                signs.append(bound.lo>0 if value>0 else bound.hi<0)
        if all(signs):box_width=width;break
    assert box_width is not None
    out=dict(controls=[str(x) for x in controls],M=[[str(x) for x in row] for row in M],T=[[str(x) for x in row] for row in T],
        determinant_coefficients_ascending=[str(x) for x in poly],root_brackets=[[str(x) for x in pair] for pair in brackets],
        endpoint_signs=[[1 if x>0 else -1 for x in pair] for pair in endpoint_values],
        velocity_lower_bound=str(velocity_lower),velocity_lower_bound_float=float(velocity_lower),
        independent_control_box_halfwidth=str(box_width),exact_inverse_root_and_box_checks=True)
    return out,T,M,poly


if __name__=='__main__':
    out,_,_,_=certificate();(Path(__file__).parent/'exact.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('T','determinant_coefficients_ascending')},indent=2))

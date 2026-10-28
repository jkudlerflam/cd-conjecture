"""Small checks of C50 lemmas. These do not compute a microscopic RG flow."""
from fractions import Fraction as F
from itertools import product
from math import exp, log
from pathlib import Path
import json
import numpy as np


def cayley_check():
    records=[]
    for r in range(2,7):
        weights=list(range(2,r+2))
        total=0
        for sequence in product(range(r),repeat=r-2):
            degrees=[1]*r
            for i in sequence: degrees[i]+=1
            term=1
            for p,d in zip(weights,degrees): term*=p**d
            total+=term
        expected=int(np.prod(weights))*sum(weights)**(r-2)
        assert total==expected
        records.append({'vertices':r,'tree_sum':total,'formula':expected})
    return records


def mul(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():
            if x&y: continue
            inversions=sum(1 for i in range(8) if x>>i&1
                             for j in range(i) if y>>j&1)
            z=x|y
            out[z]=out.get(z,F(0))+(-1)**inversions*c*d
    return {k:v for k,v in out.items() if v}


def add(a,b):
    out=dict(a)
    for k,v in b.items(): out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}


def scale(a,c): return {k:v*c for k,v in a.items() if v*c}


def grassmann_check():
    cov={(4,5):F(2,5),(4,7):F(1,10),(6,5):F(1,10),(6,7):F(3,10)}
    cov.update({(b,a):-v for (a,b),v in list(cov.items())})
    def wick(ids):
        if not ids: return F(1)
        if len(ids)%2: return F(0)
        return sum(((-1)**(j-1)*cov.get((ids[0],ids[j]),F(0))
                    *wick(ids[1:j]+ids[j+1:]) for j in range(1,len(ids))),F(0))
    def expect(poly):
        out={}
        for mask,c in poly.items():
            ids=[i for i in range(4,8) if mask>>i&1]
            q=mask&15
            out[q]=out.get(q,F(0))+c*wick(ids)
        return {k:v for k,v in out.items() if v}
    def gen(i): return {1<<i:F(1)}
    fields=[add(gen(i),gen(i+4)) for i in range(4)]
    rho1=mul(fields[0],fields[1]);rho2=mul(fields[2],fields[3])
    lam=F(1,10**12)
    potential=scale(add(mul(rho1,rho2),scale(rho1,F(1,5))),lam)
    exponential={0:F(1)};term={0:F(1)}
    for r in range(1,6):
        term=scale(mul(term,potential),F(1,r))
        exponential=add(exponential,term)
    averaged=expect(exponential)
    z=averaged.get(0,F(0)); assert z>0
    nilpotent=scale(averaged,1/z);nilpotent.pop(0)
    effective={};power={0:F(1)}
    for r in range(1,5):
        power=mul(power,nilpotent)
        effective=add(effective,scale(power,F((-1)**(r+1),r)))
    radius=8*exp(2)
    epsilon=sum(abs(float(c))*radius**mask.bit_count() for mask,c in potential.items())
    output_norm=abs(log(float(z)))+sum(abs(float(c)) for c in effective.values())
    bound=epsilon/(1-epsilon)
    assert epsilon<1 and output_norm<=bound
    # An exact six-field first-order contraction.
    six=mul(mul(mul(gen(0),gen(1)),mul(gen(2),gen(3))),mul(gen(4),gen(5)))
    contracted=expect(six)
    assert contracted=={15:F(2,5)}
    return {'input_H_norm':epsilon,'output_norm':output_norm,'bound':bound,
            'six_field_to_quartic_coefficient':str(contracted[15]),
            'exact_partition_function':str(z)}


def recurrence_check():
    A,b,eps,kappa,theta=1.7,2.0,.01,.8,.3
    C=2**(-(kappa-theta))/((1-2**(-kappa))*(1-2**(-(kappa-theta))))
    assert b*eps*C<.5
    B={}
    for h in range(0,-161,-1):
        history=sum(2**(-kappa*(k-h))*sum(B[j] for j in range(h+1,k+1))
                    for k in range(h+1,1))
        B[h]=A*eps*2**(theta*h)+b*eps*history
    ratio=max(B[h]/(2*A*eps*2**(theta*h)) for h in B)
    assert ratio<=1
    partial=sum(B.values())
    bound=2*A*eps/(1-2**(-theta))
    assert partial<=bound
    return {'condition':b*eps*C,'max_bound_ratio':ratio,
            'sum_161_increments':partial,'geometric_bound':bound}


def root_check():
    roots=np.array([-1.1,.13,1.9])
    speeds=np.array([1.,-.07,1.e-4])
    eps=1.e-8
    A=np.array([[.1,.03,-.02],[.02,-.08,.01],[.01,.02,.05]])
    c=np.array([.2,-.1,.3])
    def nu(p):return eps*(c+A@np.tanh(p-roots))
    p=roots.copy()
    for i in range(50):
        z=-nu(p)/speeds
        q=2*z/(1+np.sqrt(1+.4*z))
        pn=roots+q
        if np.max(np.abs(pn-p))<1.e-15:
            p=pn;break
        p=pn
    q=p-roots
    residual=float(np.max(np.abs(speeds*(q+.1*q*q)+nu(p))))
    shift=float(np.max(np.abs(q)))
    Cnu=float(np.max(np.abs(c)+np.sum(np.abs(A),axis=1)))
    bound=2*Cnu*eps/np.min(np.abs(speeds))
    assert residual<1.e-14 and shift<=bound
    return {'velocity_ratio':1.e4,'residual':residual,'shift':shift,'shift_bound':bound}


def main():
    result={'status':'small controls only, no all-scale RG or microscopic calculation',
            'weighted_cayley':cayley_check(),'finite_grassmann':grassmann_check(),
            'conditional_recurrence':recurrence_check(),'conditional_root_map':root_check()}
    Path('checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()

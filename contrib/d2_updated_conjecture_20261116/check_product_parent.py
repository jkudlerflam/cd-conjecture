"""Exact product-projector continuation of the deeper cancellation.

The telescoping bond gauge is checked explicitly. No c or gap is fitted.
"""
import json
from pathlib import Path
import numpy as np


def bond(r):
    n = np.diag([0.,1.]); I = np.eye(2)
    bd = np.array([[0.,0.],[1.,0.]]); b=bd.T; X=bd+b
    V=(1-r*r)**2; lam=r*(1-r*r); w=p=-r*r; mu=-2*r*r*(1-r*r)
    h = (V*np.kron(n,n)+lam*(np.kron(n,X)-np.kron(X,n))
         +w*(np.kron(b,bd)+np.kron(bd,b))
         +p*(np.kron(bd,bd)+np.kron(b,b))
         -mu/2*(np.kron(n,I)+np.kron(I,n)))
    psi = np.array([-r*r,-r,r,1.])
    gauge = r**3*(np.kron(I,X)-np.kron(X,I))
    return h,psi,gauge


def main():
    rng = np.random.default_rng(20261117)
    error=0.
    for r in rng.uniform(0.02,0.95,200):
        h,psi,gauge = bond(r)
        err=float(np.max(abs(h+r**4*np.eye(4)+gauge-np.outer(psi,psi))))
        error=max(error,err)
        assert err < 1e-12
        V=(1-r*r)**2; lam=r*(1-r*r)
        assert abs(-lam*lam/V+r*r)<1e-12
        assert abs(V-np.sqrt(V*V+4*lam*lam)+2*r*r*(1-r*r))<1e-12
    rows=[]
    r=0.3
    h,psi,_=bond(r)
    for L in range(4,11):
        H=np.zeros((1 << L,1 << L))
        for s in range(1 << L):
            for j in range(L):
                k=(j+1)%L; loc=2*((s >> j)&1)+((s >> k)&1)
                base=s & ~((1 << j)|(1 << k))
                for out in range(4):
                    H[base|((out >> 1) << j)|((out & 1) << k),s]+=h[out,loc]
        shifted=H+L*r**4*np.eye(1 << L)
        ev=np.linalg.eigvalsh(shifted)
        a=np.array([1.,r])/np.sqrt(1+r*r)
        b=np.array([1.,-r])/np.sqrt(1+r*r)
        av=np.ones(1); bv=np.ones(1)
        for _ in range(L): av=np.kron(av,a); bv=np.kron(bv,b)
        residual=max(float(np.max(abs(shifted@av))),float(np.max(abs(shifted@bv))))
        dim=int(np.count_nonzero(abs(ev)<1e-9))
        assert dim == 2 and residual < 1e-12 and ev[0]>-1e-10
        rows.append({'L':L,'ground_dimension':dim,'ground_state_residual':residual,
                     'finite_size_gap':float(ev[2])})
    out={'seed':20261117,'factorization_cases':200,'max_absolute_identity_error':error,
         'ring_checks':rows,'interpretation':'Exact low-entanglement parent; no CFT or gap fit'}
    Path(__file__).with_name('product_parent_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__': main()

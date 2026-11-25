"""Production quartet offsets on a complex parity cat, checked analytically.

This is an independent dense local-matrix product calculation for a
two-product-state superposition, without exponentially large dense vectors.
The random L12 dense-vector check in validate.py is complementary.
"""
import argparse
import numpy as np
from tenpy.networks.mps import MPS
from tenpy.networks.site import SpinHalfSite
from common import write_json
import measure_matrix


def check():
    L=24;rng=np.random.default_rng(20261126)
    v=rng.normal(size=(L,2))+1j*rng.normal(size=(L,2));v/=np.linalg.norm(v,axis=1)[:,None]
    z=np.diag([1.,-1.]);sp=np.array([[0.,1.],[0.,0.]]);sm=sp.T
    O={'Sm':sm,'Sp':sp,'Sigmaz':z}
    w=v@z;states=[v,w]
    site=SpinHalfSite(conserve=None)
    a=MPS.from_product_state([site]*L,v,dtype=complex)
    b=MPS.from_product_state([site]*L,w,dtype=complex)
    psi=a.add(b,1,1);psi.norm=1.
    def ev(ops):
        mats=[np.eye(2,dtype=complex) for _ in range(L)]
        for name,j in ops: mats[j]=mats[j]@O[name]
        return sum(np.prod([np.vdot(x[j],mats[j]@y[j]) for j in range(L)]) for x in states for y in states)
    norm=ev([])
    err=0.;count=0;by_layout={}
    for centered,refs in ((False,[0,1,2]),(True,[4,5])):
        m=measure_matrix.measure(psi,refs,rmin=10,rmax=13,progress=False,centered=centered)
        localerr=0.;localcount=0
        for i,row in m['rows'].items():
            C=np.array(row['C']);C=C[...,0]+1j*C[...,1]
            for ai,an in enumerate(m['names']):
                aa=[(n,int(i)+k) for n,k in m['operators'][an]]
                for bi,bn in enumerate(m['names']):
                    for t,r in enumerate(row['r']):
                        bb=[(n,int(i)+r+k) for n,k in m['operators'][bn]]
                        ex=ev(aa+bb)/norm-ev(aa)*ev(bb)/norm**2
                        localerr=max(localerr,abs(C[t,ai,bi]-ex));localcount+=1
        by_layout[str(centered)]={'entries':localcount,'max_abs_error':float(localerr)}
        err=max(err,localerr);count+=localcount
    assert err<1e-12,err
    return {'passed':True,'L':L,'seed':20261126,'entries':count,'max_abs_error':float(err),
            'state':'complex random product state projected onto even parity','exact_production_offsets':True,
            'by_centered_layout':by_layout}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    r=check();write_json(a.output,r);print(r)

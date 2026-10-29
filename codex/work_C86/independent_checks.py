"""Exact Fraction inverses independently check selected Jacobi enclosures."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from fractions import Fraction as Q
import numpy as np
from certify_obstructions import inverse, path_certificate, two_site_certificate


def exact(n,D=0,c=0,d=0):
    s=n+1
    H=np.diag([4]*n+[4+D]).astype(object)
    for j in range(n):H[j,j+1]=H[j+1,j]=-1
    R=inverse(np.kron(H,np.eye(s,dtype=int))+np.kron(np.eye(s,dtype=int),H))
    left=np.eye(s,dtype=complex)[0];left[1]=c
    right=np.eye(s,dtype=complex)[-1];right[-2]=d
    b=np.outer(right,left).ravel()
    xr=R@np.array([Q(float(z.real)) for z in b],dtype=object)
    xi=R@np.array([Q(float(z.imag)) for z in b],dtype=object)
    C=np.zeros((n,n),dtype=object)
    for x in [xr,xi]:
        U=np.zeros((s*s,n),dtype=object);V=U.copy()
        for j in range(n):
            U[j*s:(j+1)*s,j]=x[j*s:(j+1)*s]
            V[j::s,j]=x[j::s]
        C+=U.T@R@V
    S=(C+C.T)/2
    return S


def dec(a):return Q(int(a[0]),int(a[1]))
def enc(a):return [str(a.numerator),str(a.denominator)]


def main():
    out=[]
    for n,D in [(2,1000),(3,1)]:
        S=exact(n,D);cert=path_certificate(n,D)
        indices=[0,n-1]
        for a,i in enumerate(indices):
            for b,j in enumerate(indices):
                assert dec(cert['hermitian_minor_lower'][a][b])<=S[i,j]
                assert S[i,j]<=dec(cert['hermitian_minor_upper'][a][b])
        det=S[0,0]*S[-1,-1]-S[0,-1]**2
        upper=dec(cert['hermitian_minor_determinant_upper'])
        assert det<=upper<0
        out.append(dict(status='proved_exact_cross_check',family='point_sources',n=n,D=D,
            exact_determinant=enc(det),upper_bound=enc(upper)))
    for c,d in [(1j,(1-1j)/2),(-1,1j)]:
        n=3;S=exact(n,c=c,d=d);cert=two_site_certificate(n,c,d)
        det=S[0,0]*S[-1,-1]-S[0,-1]**2
        upper=dec(cert['determinant_upper'])
        assert det<=upper<0
        out.append(dict(status='proved_exact_cross_check',family='extended_sources',n=n,
            c=[float(c.real),float(c.imag)],d=[float(d.real),float(d.imag)],
            exact_determinant=enc(det),upper_bound=enc(upper)))
    json.dump(out,open('independent_checks.json','w'),indent=2)
    print(json.dumps(dict(status='passed_exact_inverse_cross_checks',cases=len(out),max_product_dimension=16)))


if __name__=='__main__':main()

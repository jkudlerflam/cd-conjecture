"""C84-style floor-rounded resolvents, with exact rejection witnesses.

The Jacobi enclosure is adapted from C84/certify_tree.py. This tests a
single path rather than its folded graph. No sibling module is imported.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import json
from fractions import Fraction as Q
import numpy as np


def encode(q):
    q = Q(q)
    return [str(q.numerator), str(q.denominator)]


def path_certificate(n, end_mass, bits=256):
    s = n+1
    scale = 1 << bits
    diag = np.array([4]*n+[4+end_mass], dtype=object)
    den = diag[:, None]+diag[None, :]
    steps = bits+24
    def solve(B):
        X = np.zeros(B.shape, dtype=object)
        div = den if B.ndim == 2 else den[:, :, None]
        for _ in range(steps):
            Y = B.copy()
            Y[1:] += X[:-1]
            Y[:-1] += X[1:]
            Y[:, 1:] += X[:, :-1]
            Y[:, :-1] += X[:, 1:]
            X = Y//div
        return X
    b = np.zeros((s, s), dtype=object)
    b[-1, 0] = scale
    x = solve(b)
    V = np.zeros((s, s, n), dtype=object)
    for j in range(n):
        V[:, j, j] = x[:, j]
    y = solve(V)
    C = np.zeros((n, n), dtype=object)
    Cu = C.copy()
    for i in range(n):
        for j in range(n):
            C[i, j] = sum(x[i, k]*y[i, k, j] for k in range(s))
            Cu[i, j] = sum((x[i, k]+3)*(y[i, k, j]+4) for k in range(s))
    L, U = C+C.T, Cu+Cu.T
    i, j = 0, n-1
    det_upper_num = U[i, i]*U[j, j]-L[i, j]**2
    assert L[i, i]>0 and L[j, j]>0 and det_upper_num<0
    determinant_upper=Q(int(det_upper_num),4*scale**4)
    norm_upper=Q(int(max(U[i,i]+U[i,j],U[j,j]+U[i,j])),2*scale**2)
    result = dict(status='proved_integer_interval_obstruction', n=n,
        scalar_vertices=s, product_dimension=s*s, end_mass=end_mass,
        bits=bits, iterations=steps, active_sites=list(range(n)),
        principal_indices=[i, j],
        hermitian_minor_determinant_upper=encode(determinant_upper),
        positive_and_negative_sign_gap_lower=encode(-determinant_upper/norm_upper),
        hermitian_minor_lower=[[encode(Q(int(L[r,c]),2*scale**2)) for c in [i,j]] for r in [i,j]],
        hermitian_minor_upper=[[encode(Q(int(U[r,c]),2*scale**2)) for c in [i,j]] for r in [i,j]],
        outcome='indefinite_hopping_rejected_before_root_stage')
    return result


def inverse(A):
    n=len(A)
    B=[[Q(v) for v in row]+[Q(i==j) for j in range(n)] for i,row in enumerate(A)]
    for i in range(n):
        k=next(k for k in range(i,n) if B[k][i])
        B[i],B[k]=B[k],B[i]
        d=B[i][i];B[i]=[x/d for x in B[i]]
        for k in range(n):
            if k==i:continue
            d=B[k][i];B[k]=[x-d*y for x,y in zip(B[k],B[i])]
    return np.array([row[n:] for row in B],dtype=object)


def two_vertex_checks():
    records=[]
    for a in [Q(3,2),Q(2),Q(4),Q(9)]:
        for b in [Q(1),Q(3),Q(5)]:
            H=np.array([[a,-1],[-1,b]],dtype=object)
            A=np.kron(H,np.eye(2,dtype=int))+np.kron(np.eye(2,dtype=int),H)
            R=inverse(A);x=R[:,2].reshape(2,2)
            U=np.zeros((4,2),dtype=object);V=U.copy()
            for j in range(2):
                U[2*j:2*j+2,j]=x[j,:];V[j::2,j]=x[:,j]
            C=U.T@R@V;S=(C+C.T)/2
            det=S[0,0]*S[1,1]-S[0,1]**2
            z=a*b;s=a*a+b*b
            P=s*s+(8*z*z+z-1)*s+16*z**4-16*z**3+18*z*z-6*z
            formula=-P/(64*(a+b)**6*(z-1)**4)
            assert det==formula and det<0
            asym=C[1,0]-C[0,1]
            asym_formula=(s+2*z*(2*z-1))/(4*(a+b)**3*(z-1)**2)
            assert asym==asym_formula and asym>0
            records.append(dict(a=encode(a),b=encode(b),determinant=encode(det),asymmetry=encode(asym)))
    return records


def two_site_certificate(n, c, d, bits=256):
    s=n+1;scale=1<<bits;steps=bits+24
    def solve(B):
        X=np.zeros(B.shape,dtype=object)
        for _ in range(steps):
            Y=B.copy();Y[1:]+=X[:-1];Y[:-1]+=X[1:]
            Y[:,1:]+=X[:,:-1];Y[:,:-1]+=X[:,1:]
            X=Y//8
        return X
    left=np.eye(s,dtype=complex)[0];left[1]=c
    right=np.eye(s,dtype=complex)[-1];right[-2]=d
    source=np.outer(right,left)
    br=np.array([[int(Q(float(z.real))*scale) for z in row] for row in source],dtype=object)
    bi=np.array([[int(Q(float(z.imag))*scale) for z in row] for row in source],dtype=object)
    xr,xi=solve(br),solve(bi)
    vr=np.zeros((s,s,n),dtype=object);vi=vr.copy()
    for j in range(n):
        vr[:,j,j]=xr[:,j];vi[:,j,j]=xi[:,j]
    yr,yi=solve(vr),solve(vi)
    C=np.zeros((n,n),dtype=object);E=C.copy()
    for i in range(n):
        for j in range(n):
            C[i,j]=sum(xr[i,k]*yr[i,k,j]+xi[i,k]*yi[i,k,j] for k in range(s))
            E[i,j]=sum(4*(abs(xr[i,k])+abs(xi[i,k]))+
                3*(abs(yr[i,k,j])+abs(yi[i,k,j]))+24 for k in range(s))
    S=C+C.T;E=E+E.T
    i,j=0,n-1
    diagonal_product_upper=max((S[i,i]+si*E[i,i])*(S[j,j]+sj*E[j,j])
        for si in [-1,1] for sj in [-1,1])
    off_diagonal_abs_lower=max(abs(S[i,j])-E[i,j],0)
    upper=diagonal_product_upper-off_diagonal_abs_lower**2
    assert upper<0
    determinant_upper=Q(int(upper),4*scale**4)
    norm_upper=Q(int(max(abs(S[i,i])+E[i,i]+abs(S[i,j])+E[i,j],
        abs(S[j,j])+E[j,j]+abs(S[i,j])+E[i,j])),2*scale**2)
    return dict(status='proved_integer_interval_obstruction',n=n,vertices=s,
        source_left=[encode(Q(float(c.real))),encode(Q(float(c.imag)))],
        source_right=[encode(Q(float(d.real))),encode(Q(float(d.imag)))],
        bits=bits,iterations=steps,principal_indices=[i,j],
        determinant_upper=encode(determinant_upper),
        positive_and_negative_sign_gap_lower=encode(-determinant_upper/norm_upper),
        outcome='no_definite_hopping_for_any_common_port_phase')


def main():
    amps=[0,1,-1,1j,-1j,(1+1j)/2,(1-1j)/2]
    out=dict(status='proved_exact_finite_checks',
        two_vertex_formula_checks=two_vertex_checks(),
        single_rail_certificates=[path_certificate(n,D) for n in [2,3,6,9,10] for D in [0,1,10,1000]],
        two_site_source_certificates=[two_site_certificate(n,c,d) for n in [9,10] for c in amps for d in amps])
    json.dump(out,open('obstruction_certificate.json','w'),indent=2)
    print(json.dumps(dict(status=out['status'],two_vertex_cases=len(out['two_vertex_formula_checks']),
        single_rail_cases=len(out['single_rail_certificates']),
        two_site_source_cases=len(out['two_site_source_certificates']),max_product_dimension=121)))


if __name__=='__main__':
    main()

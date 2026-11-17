"""Exact rational/polynomial identities; requires sympy as well as numpy/scipy.

The literal spin matrices enter before any claimed effective formula.
This is finite-ring verification, not a thermodynamic convergence theorem.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import sympy as sp
from scipy import sparse

from algebra import bond, chain, occupations, projectors


def exact(A):
    A=A.tocoo()
    def rat(x):
        return sp.Rational(float(x)).limit_denominator(1000000)
    return sp.MutableSparseMatrix(A.shape[0], A.shape[1],
                                 {(int(i),int(j)):rat(z.real)+sp.I*rat(z.imag)
                                  for i,j,z in zip(A.row,A.col,A.data)})


def take(A, rows, cols=None):
    return A.extract(list(map(int,rows)), list(map(int,rows if cols is None else cols)))


def expand(A):
    return A.applyfunc(sp.expand)


def zero(A):
    return all(sp.expand(x)==0 for x in A.todok().values())


def fourth(L, dense):
    x,y,r,s=sp.symbols('x y r s', real=True)
    q=x+sp.I*y; f=r+sp.I*s
    ns=occupations(L); E=np.sum(ns*np.roll(ns,-1,axis=1),axis=1)
    D2=np.sum(ns*np.roll(ns,-2,axis=1),axis=1)
    P=np.flatnonzero(E==0); Q=np.flatnonzero(E>0)
    constraint=L-3*np.sum(ns,axis=1)+D2 if dense else D2
    keep=np.flatnonzero(constraint[P]==0); remove=np.flatnonzero(constraint[P]>0)
    F=exact(chain(L,bond(a=1,b=-1)))
    W=(exact(chain(L,bond(u=-4 if dense else 2,t=-1)))
       +x*exact(chain(L,bond(p=1)))+y*exact(chain(L,bond(p=1j))))
    V3=r*exact(chain(L,bond(f=1)))+s*exact(chain(L,bond(f=1j))) if dense else sp.zeros(1<<L)
    if dense: W+=2*L*sp.eye(1<<L)
    R=sp.diag(*[sp.Rational(1,int(e)) for e in E[Q]])
    A=take(F,Q,P); B=take(W,Q,P); C=take(F,Q); D=take(W,Q)
    RA=R*A; RB=R*B; CRA=C*RA
    S2=take(W,P)-A.adjoint()*RA
    S3=expand(take(V3,P)+RA.adjoint()*C*RA-A.adjoint()*RB-B.adjoint()*RA)
    S4=expand(-B.adjoint()*RB+CRA.adjoint()*RB+RB.adjoint()*CRA
              +RA.adjoint()*D*RA-CRA.adjoint()*R*CRA
              -RA.adjoint()*take(V3,Q,P)-take(V3,P,Q)*RA)
    R2=sp.diag(*[sp.Rational(1,2*int(e)) for e in constraint[P[remove]]])
    K=expand(take(S4,keep)-take(S3,keep,remove)*R2*take(S3,remove,keep))
    target=P[keep]; nn=ns[target]; N=np.sum(nn,axis=1)
    D3=np.sum(nn*np.roll(nn,-3,axis=1),axis=1)
    D4=np.sum(nn*np.roll(nn,-4,axis=1),axis=1)
    T=take(exact(chain(L,bond(t=1))),target)
    Ti=take(exact(chain(L,bond(t=1j))),target)
    qq=x*x+y*y; ff=r*r+s*s
    if dense:
        c=sp.Rational(1,2)-x+(qq+ff)/6
        u=-sp.Rational(7,2)+3*x-5*qq/6-2*ff/3
        a=((x-1)**2+y*y)/2-ff/6
        tre=((x-1)**2+y*y-ff)/4
        tim=((3+x)*s-y*r)/2
        expected=sp.diag(*[c*L+u*int(n)+a*int(d4) for n,d4 in zip(N,D4)])+tre*T+tim*Ti
    else:
        expected=sp.diag(*[-L*qq+(4+4*x+2*qq)*int(n)
                          +sp.Rational(2,3)*((x-3)**2+y*y)*int(d3)
                          +((x-1)**2+y*y)*int(d4)
                          for n,d3,d4 in zip(N,D3,D4)])-2*(1+x)*T
    assert zero(S2-sp.diag(*[2*int(v) for v in constraint[P]]))
    assert zero(take(S3,keep))
    assert zero(K-expected)
    return {'L':L,'branch':'dense' if dense else 'dilute',
            'dimension':len(target),'free_real_symbols':4 if dense else 2,
            'every_matrix_entry_exact':True}


def schur(coeffs,keep):
    remove=sorted(set(range(coeffs[0].rows))-set(keep))
    R=sp.diag(*[1/coeffs[0][j,j] for j in remove])
    B=[take(c,remove,keep) for c in coeffs]; C=[take(c,remove) for c in coeffs]
    X=[sp.zeros(len(remove),len(keep))]; out=[take(coeffs[0],keep)]
    for n in range(1,len(coeffs)):
        v=B[n]
        for k in range(1,n):v=v-C[k]*X[n-k]
        X.append(expand(R*v))
        v=take(coeffs[n],keep)
        for k in range(1,n):v=v-B[k].adjoint()*X[n-k]
        out.append(expand(v))
    return out


def dense_sixth_exact(L,which):
    """Exact full matrices at the isolated point and circle endpoints."""
    if which=='generic':
        x,y,r,s=sp.symbols('x y r s',real=True);q=x+sp.I*y;f=r+sp.I*s
    elif which=='isolated': q=sp.Integer(1); f=sp.Integer(0)
    elif which=='circle_zero': q=sp.Integer(0); f=sp.sqrt(3)
    else: q=sp.Integer(3); f=2*sp.I*sp.sqrt(3)
    qq=sp.expand(q*sp.conjugate(q)); ff=sp.expand(f*sp.conjugate(f)); x=sp.re(q)
    c=sp.Rational(1,2)-x+(qq+ff)/6
    u=-sp.Rational(7,2)+3*x-5*qq/6-2*ff/3
    tre=((x-1)**2+sp.im(q)**2-ff)/4
    tim=sp.im((3+sp.conjugate(q))*f)/2
    dim=1<<L; identity=sp.eye(dim)
    H=[exact(chain(L,bond(v=1))),exact(chain(L,bond(a=1,b=-1))),
       exact(chain(L,bond(u=-4,t=-1)))+sp.re(q)*exact(chain(L,bond(p=1)))
       +sp.im(q)*exact(chain(L,bond(p=1j)))+2*L*identity,
       sp.re(f)*exact(chain(L,bond(f=1)))+sp.im(f)*exact(chain(L,bond(f=1j))),
       -u*exact(chain(L,bond(u=1)))-tre*exact(chain(L,bond(t=1)))
       -tim*exact(chain(L,bond(t=1j)))-c*L*identity,
       sp.zeros(dim),sp.zeros(dim)]
    P=list(map(int,projectors(L))); S=schur(H,P)
    keep=[i for i in range(len(P)) if S[2][i,i]==0]
    K=schur(S[2:],keep); words=[P[i] for i in keep]
    if which!='generic': assert all(zero(m) for m in K[:4])
    source=sum(1<<j for j in (0,2,4,7)); dest=sum(1<<j for j in (0,3,5,7))
    for j in range(10,L-1,2):source+=1<<j;dest+=1<<j
    actual=sp.simplify(K[4][words.index(dest),words.index(source)])
    if which=='generic':
        before=(1<<0)+(1<<3)+(1<<6)
        after=(1<<0)+(1<<2)+(1<<4)+(1<<6)
        gamma=sp.expand(K[3][words.index(after),words.index(before)])
        expected_gamma=(q+3)*sp.conjugate(f)/6+(sp.conjugate(q)-1)*f/2
        cone=sp.groebner([r*r+s*s-3*((x-1)**2+y*y)],r,s,x,y,domain=sp.QQ)
        assert cone.reduce(sp.expand(sp.re(gamma-expected_gamma)))[1]==0
        assert cone.reduce(sp.expand(sp.im(gamma-expected_gamma)))[1]==0
        from check import conversion
        raising=exact(sparse.csr_matrix(conversion(L,np.array(words),1)))
        phase=exact(sparse.csr_matrix(conversion(L,np.array(words),1j)))
        difference=expand(K[3]-sp.re(expected_gamma).expand(complex=True)*raising
                          -sp.im(expected_gamma).expand(complex=True)*phase)
        for value in difference.todok().values():
            assert cone.reduce(sp.expand(sp.re(value)))[1]==0
            assert cone.reduce(sp.expand(sp.im(value)))[1]==0
        ideal=sp.groebner([y*y+x*x-3*x,r*r-(3-x)*(x+1),s*s-x*(x+1),
                           x*r+y*s,(3-x)*s+y*r],r,s,y,x,domain=sp.QQ)
        expected=(15*x*x+4*x-3)/16+sp.I*(5*x-1)*s/4
        assert ideal.reduce(sp.expand(sp.re(actual-expected)))[1]==0
        assert ideal.reduce(sp.expand(sp.im(actual-expected)))[1]==0
        assert sp.simplify(actual.subs({x:1,y:0,r:0,s:0})+sp.Rational(1,2))==0
        return {'L':L,'case':'generic','fifth_conversion':str(gamma),
                'sixth_witness_unrestricted':str(sp.expand(actual)),
                'circle_reduced_witness':str(expected),'exact_polynomial_reductions_pass':True,
                'entire_fifth_operator_checked':True}
    expected=-sp.Rational(1,2) if which=='isolated' else ((15*x*x+4*x-3)/16+sp.I*(5*x-1)*sp.im(f)/4)
    assert sp.simplify(actual-expected)==0
    return {'L':L,'case':which,'sixth_witness':str(actual),'all_lower_coefficients_exactly_zero':True}


def dilute_sixth_exact(L,q):
    from check_dilute_sixth import local_sixth
    dim=1<<L
    H=[exact(chain(L,bond(v=1))),exact(chain(L,bond(a=1,b=-1))),
       exact(chain(L,bond(u=2,t=-1,p=q))),sp.zeros(dim),
       exact(chain(L,bond(u=-(4+4*q+2*q*q),t=2+2*q)))+L*q*q*sp.eye(dim),
       sp.zeros(dim),sp.zeros(dim)]
    P=list(map(int,projectors(L)));S=schur(H,P)
    keep2=[i for i in range(len(P)) if S[2][i,i]==0]
    K=schur(S[2:],keep2)
    keep4=[i for i in range(len(keep2)) if K[2][i,i]==0]
    last=schur(K[2:],keep4);words=[P[keep2[j]] for j in keep4]
    assert zero(last[0]) and zero(last[1])
    expected=exact(sparse.csr_matrix(local_sixth(L,np.array(words),q)))
    assert zero(last[2]-expected)
    return {'L':L,'q':q,'dimension':len(words),'sixth_full_matrix_exact':True}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--sixth',action='store_true')
    parser.add_argument('--generic-sixth',action='store_true')
    parser.add_argument('--dilute-sixth',action='store_true')
    args=parser.parse_args();start=time.perf_counter()
    rows=[]
    for dense in (False,True):
        rows.append(fourth(9,dense)); print(rows[-1],flush=True)
    six=[]
    if args.sixth:
        for which in ('isolated','circle_zero','circle_three'):
            six.append(dense_sixth_exact(9,which));print(six[-1],flush=True)
    if args.generic_sixth:
        six.append(dense_sixth_exact(9,'generic'));print(six[-1],flush=True)
    dilute=[]
    if args.dilute_sixth:
        for q in (0,-1,1,3):
            dilute.append(dilute_sixth_exact(10,q));print(dilute[-1],flush=True)
    out={'sympy_version':sp.__version__,'fourth_order':rows,'sixth_order':six,
         'dilute_sixth_order':dilute,
         'runtime_seconds':time.perf_counter()-start}
    Path(__file__).with_name('symbolic_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()

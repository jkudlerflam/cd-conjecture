"""Bounded one-worker checks of the actual loop and analytic coefficients."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from itertools import product,combinations_with_replacement
import json,time
import numpy as np
import scipy.linalg as la
import mpmath as mp
from loop import carrier,coefficients,low_limit,roots_and_response
from moments import criterion

HERE=Path(__file__).resolve().parent
mp.mp.dps=75


def clock(n,q):
    p=n+1
    H=np.diag([q*q]+[1+q*q]*p)-q*(np.eye(p+1,k=1)+np.eye(p+1,k=-1))
    e,v=la.eigh(H)
    return 4+q*q-e[0],np.array([sum(v[j+1:,0]**2)/q**2 for j in range(n)])


def parity():
    n=2;D=20.;w=.5;eta=1e-4;t=eta/D**(n+1)
    T0,M1=coefficients(n,D,w)
    T,M=coefficients(n,D,w,td=t);Tm,Mm=coefficients(n,D,w,td=-t)
    assert np.array_equal(T,Tm) and np.array_equal(M,-Mm)
    error=la.norm(M/t-M1)/la.norm(M1)
    assert error<1e-6
    _,single=coefficients(1,D,w)
    assert la.norm(single)<1e-20
    return dict(n=n,D=D,eta=eta,onsite_derivative_relative_error=error,
                hopping_even=True,onsite_odd=True,n1_derivative_norm=float(la.norm(single)))


def finite_delay():
    reports=[]
    for n in (2,3):
        p=n+1;w=.5;q=.0625;a,gamma=clock(n,q)
        G,B,_,Q=low_limit(n,w,a)
        F0=np.diag(gamma);nu=.1/max(la.eigvalsh(F0,G));F=nu*F0
        e,X=la.eigh(F,G);e=e[::-1];X=X[:,::-1]
        QQ=X.T@Q@X;BB=X.T@B@X
        expected=np.array([-16*(p+1)*sum(QQ[j,l]*BB[j,l]/(e[j]-e[l])
                                         for l in range(n) if l!=j) for j in range(n)])
        rows=[]
        for D in (1e4,1e5,1e6):
            T,M1=coefficients(n,D,w,a);scale=D**(2*p+2)
            roots=roots_and_response(T*scale,M1*scale/D**p,F)
            response=np.array([z[1]*D for z in roots])
            error=float(la.norm(response-expected)/la.norm(expected))
            qerr=float(la.norm(scale*D*(T-T.T)/(4*(p+1))-Q)/la.norm(Q))
            berr=float(la.norm(M1*scale/D**p-2j*B)/la.norm(B))
            rows.append(dict(D=D,D_times_pair_response=response.tolist(),
                             response_relative_error=error,delay_relative_error=qerr,
                             loop_source_relative_error=berr))
        assert rows[-1]['response_relative_error']<1e-3
        assert rows[-1]['delay_relative_error']<1e-3
        assert rows[-1]['loop_source_relative_error']<1e-3
        assert rows[-1]['response_relative_error']<rows[0]['response_relative_error']/20
        reports.append(dict(n=n,w=w,q=q,a=a,nu=nu,expected=expected.tolist(),checks=rows))
    return reports


def shortest_walks():
    reports=[]
    for n in (2,3):
        for w in (.02,.01):
            G,B,_,Q=low_limit(n,w,4.)
            expected=(np.triu(np.ones((n,n)),1)-np.tril(np.ones((n,n)),-1))/8**(n+3)
            error=float(la.norm(B/w**n-expected)/la.norm(expected))
            assert error<1e-3
            reports.append(dict(n=n,w=w,B_shortest_walk_relative_error=error))
    return reports


def exact_moments():
    return [dict(n=n,L=list(map(str,criterion(n)[-1])),
                 exact=True,full_distinctness_criterion=all(x!=0 for x in criterion(n)[-1][1:]))
            for n in (2,3,4)]


def moment_asymptotics():
    reports=[]
    convert=lambda M:mp.matrix([[mp.mpf(x.numerator)/x.denominator for x in row] for row in M])
    for n in (2,3,4):
        C,Q,J,_,_,L=criterion(n);C,Q,J=map(convert,(C,Q,J))
        invL=mp.cholesky(C)**-1
        expected=[mp.mpf(x.numerator)/x.denominator for x in L]
        expected.append(-expected[-1])
        errors=[]
        for u in (mp.mpf(1)/16,mp.mpf(1)/64):
            E=mp.diag([u**(4*j) for j in range(n)])
            e,U=mp.eigsy(invL*E*invL.T)
            order=list(range(n-1,-1,-1))
            e=[e[j] for j in order]
            X=invL.T*mp.matrix([[U[i,j] for j in order] for i in range(n)])
            QQ=X.T*Q*X;BB=X.T*J*X
            values=[sum(QQ[j,l]*BB[j,l]/(e[j]-e[l]) for l in range(n) if l!=j)
                    *u**(4*min(j,n-2)) for j in range(n)]
            error=max(abs(values[j]-expected[j]) for j in range(n))
            errors.append(error)
        assert errors[1]<errors[0]/100
        reports.append(dict(n=n,limiting_coefficients=list(map(str,expected)),
                            u_values=['1/16','1/64'],maximum_absolute_errors=list(map(str,errors))))
    return reports


def mp_coefficients(n,D,w,a,t=mp.mpf(0)):
    p=n+1;s=2*p+1
    H=mp.diag([a]*p+[a+D]*(p+1));J=mp.matrix(s)
    for j in range(p-1):H[j,j+1]=H[j+1,j]=-w
    for j in range(p,2*p):H[j,j+1]=H[j+1,j]=-1
    H[0,p]=H[p,0]=-1;J[-1,p-1]=-1;J[p-1,-1]=1
    H+=1j*t*J
    L=mp.matrix(s*s);JJ=mp.matrix(s*s)
    for i,j,k,l in product(range(s),repeat=4):
        L[i*s+j,k*s+l]=H[i,k]*(j==l)+(i==k)*H[j,l]
        JJ[i*s+j,k*s+l]=J[i,k]*(j==l)+(i==k)*J[j,l]
    R=L**-1;x=R[:,(s-1)*s]
    def restrict(x):
        U=mp.matrix(s*s,n);V=mp.matrix(s*s,n)
        for j in range(n):
            for b in range(s):U[j*s+b,j]=x[j*s+b];V[b*s+j,j]=x[b*s+j]
        return U,V
    U,V=restrict(x);RU,RV=R*U,R*V
    C=U.H*RV;A=U.H*RU+V.H*RV
    T=mp.matrix([[-2*mp.re(C[i,j]) for j in range(n)] for i in range(n)])
    M=mp.matrix([[-2j*mp.im(A[i,j]) for j in range(n)] for i in range(n)])
    if t:return T,M
    YU,YV=restrict(-R*JJ*x)
    B=U.T*R*YU-YU.T*RU+V.T*R*YV-YV.T*RV-RU.T*JJ*RU-RV.T*JJ*RV
    return T,-2j*B


def high_precision_response():
    n=2;p=3;D=mp.mpf(10000);w=mp.mpf(1)/2;q=mp.mpf(1)/16;nu=mp.mpf(1)/10**6
    C=mp.diag([q*q]+[1+q*q]*p)
    for j in range(p):C[j,j+1]=C[j+1,j]=-q
    ev,v=mp.eigsy(C);a=4+q*q-ev[0]
    gamma=[sum(v[j,0]**2 for j in range(i+1,p+1))/q**2 for i in range(n)]
    F=mp.diag([nu*x for x in gamma]);scale=D**(2*p+2)
    T,M1=mp_coefficients(n,D,w,a);T*=scale;M1*=scale/D**p
    Tf=np.array(T.tolist(),complex).real;Ff=np.array(F.tolist(),complex).real
    seeds=[x[0] for x in roots_and_response(Tf,np.array(M1.tolist(),complex),Ff)]
    def roots(T,M,seeds):
        def determinant(k):return mp.re(mp.det(T*mp.exp(1j*k)+T.T*mp.exp(-1j*k)+M+F))
        return [mp.findroot(determinant,(mp.mpf(k)-mp.mpf('1e-5'),mp.mpf(k)+mp.mpf('1e-5')),
                            tol=mp.mpf('1e-65')) for k in seeds]
    k0=roots(T,mp.matrix(n),seeds)
    predicted=[]
    for k in k0:
        H=T*mp.exp(1j*k)+T.T*mp.exp(-1j*k)+F
        e,U=mp.eighe(H);j=min(range(n),key=lambda j:abs(e[j]));z=U[:,j]
        speed=mp.re((z.H*(1j*T*mp.exp(1j*k)-1j*T.T*mp.exp(-1j*k))*z)[0])
        predicted.append(-2*mp.re((z.H*M1*z)[0])/speed)
    reports=[]
    for eta in (mp.mpf('1e-6'),mp.mpf('5e-7')):
        TT,MM=mp_coefficients(n,D,w,a,eta/D**p);TT*=scale;MM*=scale
        rr=roots(TT,MM,list(k0)+[-k for k in k0])
        sums=[rr[j]+rr[n+j] for j in range(n)]
        error=max(abs(sums[j]/eta-predicted[j]) for j in range(n))
        wrap=lambda x:abs((x+mp.pi)%(2*mp.pi)-mp.pi)
        d2=min(wrap(x-y) for i,x in enumerate(rr) for y in rr[i+1:])
        pairs=list(combinations_with_replacement(range(2*n),2))
        d4=min(wrap(rr[a]+rr[b]-rr[c]-rr[d]) for i,(a,b) in enumerate(pairs) for c,d in pairs[i+1:])
        assert d4>0 and error<mp.mpf('1e-4')
        reports.append(dict(eta=str(eta),pair_sums=list(map(str,sums)),
                            pair_response_error=str(error),delta2=str(d2),delta4=str(d4)))
    assert mp.mpf(reports[1]['pair_response_error'])<mp.mpf(reports[0]['pair_response_error'])/3
    assert abs(sum(predicted))<mp.mpf('1e-50')
    return dict(n=n,decimal_precision=mp.mp.dps,D=str(D),w=str(w),q=str(q),nu=str(nu),
                exact_derivative_formula_values=list(map(str,predicted)),checks=reports,
                status='high-precision validation, not interval certification')


def main():
    start=time.monotonic();out=dict(status='bounded checks; no all-n or CFT extrapolation',workers=1)
    for name,run in [('parity',parity),('finite_delay',finite_delay),('shortest_walks',shortest_walks),
                     ('exact_moments',exact_moments),('moment_asymptotics',moment_asymptotics),
                     ('high_precision_response',high_precision_response)]:
        out[name]=run();print(name,'passed',flush=True)
        (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    out['runtime_seconds']=time.monotonic()-start
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('all checks passed in',round(out['runtime_seconds'],2),'seconds')


if __name__=='__main__':main()

"""Independent full-sector finite-ring checks. No central-charge estimates."""
import json
import time
from pathlib import Path

import numpy as np
import scipy
from scipy import sparse

from algebra import (bond, chain, dense_coefficients, dense_higher, maxabs,
                     nested_fourth, occupations, projectors, second_local, sub)


def first_kernel(L):
    names=['v','u','t','p','f','a','b']
    hs=[]
    for name in names:
        hs.append(chain(L,bond(**{name:1})))
        if name not in ('v','u'):hs.append(chain(L,bond(**{name:1j})))
    P=projectors(L); cols=[]
    for H in hs:
        K=sub(H,P).toarray();K-=K[0,0]*np.eye(len(P))
        cols.append(K.ravel())
    M=np.array(cols).T;M=np.concatenate([M.real,M.imag])
    sv=np.linalg.svd(M,compute_uv=False);rank=int(sum(sv>1e-9))
    assert rank==5
    return {'L':L,'rank':rank,'kernel_dimension':12-rank}


def fourth_local(L,q,dense=False,f=0j):
    ns=occupations(L); P=projectors(L)
    D2=np.sum(ns*np.roll(ns,-2,axis=1),axis=1)
    constraint=L-3*np.sum(ns,axis=1)+D2 if dense else D2
    P=P[constraint[P]==0]; nn=ns[P]; N=np.sum(nn,axis=1)
    D3=np.sum(nn*np.roll(nn,-3,axis=1),axis=1)
    D4=np.sum(nn*np.roll(nn,-4,axis=1),axis=1)
    if dense:
        c,u,a,t=dense_coefficients(q,f)
        out=sub(chain(L,bond(u=u,t=t)),P).toarray()+np.diag(c*L+a*D4)
    else:
        out=sub(chain(L,bond(u=4+4*q.real+2*abs(q)**2,t=-2-2*q.real)),P).toarray()
        out+=np.diag(-L*abs(q)**2+2/3*abs(q-3)**2*D3+abs(q-1)**2*D4)
    return P,out


def conversion(L,P,gamma):
    """Local 33 -> 222 tile conversion, with all spectator particles retained."""
    ns=occupations(L);pos={int(s):i for i,s in enumerate(P)}
    M=np.zeros((len(P),len(P)),complex)
    for s in P:
        for j in range(L):
            js=[(j+r)%L for r in range(7)]
            if list(ns[s,js])==[1,0,0,1,0,0,1]:
                dest=s^(1<<js[2])^(1<<js[3])^(1<<js[4])
                if dest in pos:
                    M[pos[dest],pos[s]]+=gamma;M[pos[s],pos[dest]]+=np.conj(gamma)
    return M


def eigen_check(L,q,f,order):
    P,K,H=dense_higher(L,q,f,order=6)
    lead=K[order-2].toarray(); rows=[]
    for eps in (.10,.075,.05):
        full=sum(eps**k*h for k,h in enumerate(H)).toarray()
        ev,U=np.linalg.eigh(full);m=len(P)
        left,sv,right=np.linalg.svd(U[P,:m]);polar=left@right
        effective=(polar*ev[:m])@polar.conj().T
        error=maxabs(effective/eps**order-lead)
        rows.append({'epsilon':eps,'leading_scaled_error':error,
                     'minimum_projection_overlap':float(min(sv)),
                     'band_gap_over_epsilon2':float((ev[m]-ev[m-1])/eps**2)})
        assert min(sv)>.7 and ev[m]-ev[m-1]>eps**2
    # Fifth-order leading dynamics has an O(epsilon) scaled correction;
    # the displayed real sixth-order example has an O(epsilon^2) correction.
    ratio_limit=.6 if order==5 else .4
    assert rows[-1]['leading_scaled_error']<rows[0]['leading_scaled_error']*ratio_limit
    return {'L':L,'q':[complex(q).real,complex(q).imag],
            'f':[complex(f).real,complex(f).imag],'order':order,'rows':rows}


def parent_checks(rng):
    n=np.diag([0.,1.]);I=np.eye(2);errors=[]
    for _ in range(80):
        a=.2*(rng.normal()+1j*rng.normal());c=.04*(rng.normal()+1j*rng.normal());b=1.
        psi=np.array([c,-a,a,b]);g=b*np.conj(a)+a*np.conj(c)
        h=bond(v=b*b+abs(c)**2-2*abs(a)**2,u=2*(abs(a)**2-abs(c)**2),
               t=-abs(a)**2,p=b*np.conj(c),a=g,b=-g)
        gauge=np.array([[0,c*np.conj(a)],[a*np.conj(c),0]])
        err=maxabs(h+abs(c)**2*np.eye(4)+np.kron(gauge,I)-np.kron(I,gauge)-np.outer(psi,psi.conj()))
        errors.append(err);assert err<1e-12
    return {'cases':len(errors),'max_factorization_error':max(errors)}


def main():
    start=time.perf_counter();rng=np.random.default_rng(20261116)
    first=[first_kernel(L) for L in (5,6,7)]
    seconds=[]
    for L in (7,8,9,10,11,12):
        ns=occupations(L);E=np.sum(ns*np.roll(ns,-1,axis=1),axis=1)
        P=np.flatnonzero(E==0);Q=np.flatnonzero(E>0)
        for _ in range(3):
            a,b,p,t,f=rng.normal(size=5)+1j*rng.normal(size=5);u=rng.normal()
            V1=chain(L,bond(a=a,b=b,p=p,v=.7))
            V2=chain(L,bond(u=u,t=t,f=f,a=.4+.2j,b=-.3+.6j,p=.2j,v=-.5))
            A=sub(V1,Q,P)
            K=(sub(V2,P)-A.conj().T@sparse.diags(1/E[Q])@A).toarray()
            err=maxabs(K-second_local(L,a,b,p,u,t,f));assert err<2e-11
            seconds.append({'L':L,'error':err})
    fourth=[]
    for L in (9,10,11,12,13):
        for dense in (False,True):
            for _ in range(3):
                q=rng.normal()+1j*rng.normal();f=rng.normal()+1j*rng.normal() if dense else 0j
                P,K,errors=nested_fourth(L,q,dense,source3=f,pair3=.3+.2j)
                _,expected=fourth_local(L,q,dense,f)
                err=maxabs(K-expected);assert err<5e-10
                assert errors['second_order_error']<1e-10 and errors['third_order_P2_norm']<1e-10
                fourth.append({'L':L,'dense':dense,'error':err,**errors})
    fifth=[]
    for L in (9,10,11,12):
        for _ in range(4):
            q=rng.normal()+1j*rng.normal()
            f=np.sqrt(3)*abs(q-1)*np.exp(1j*rng.uniform(0,2*np.pi))
            P,K,_=dense_higher(L,q,f,order=5)
            gamma=(q+3)*np.conj(f)/6+(np.conj(q)-1)*f/2
            err=maxabs(K[3]-conversion(L,P,gamma));assert err<1e-9
            fifth.append({'L':L,'error':err})
    sixth=[]
    for L in (9,10,12,13):
        samples=[(1.,0.)]
        for theta in np.linspace(0,np.pi,9):
            q=1.5+1.5*np.exp(1j*theta)
            f=np.sqrt(3)*abs(q-1)*np.exp(.5j*np.angle(-(q+3)/(3*(np.conj(q)-1))))
            samples.append((q,f))
        for q,f in samples:
            P,K,_=dense_higher(L,q,f);pos={int(s):i for i,s in enumerate(P)}
            source=sum(1<<j for j in (0,2,4,7));dest=sum(1<<j for j in (0,3,5,7))
            for j in range(10,L-1,2):source+=1<<j;dest+=1<<j
            actual=complex(K[4][pos[dest],pos[source]])
            expected=-.5 if f==0 else (15*q.real**2+4*q.real-3)/16+.25j*(5*q.real-1)*f.imag
            err=abs(actual-expected)
            assert maxabs(K[2])<1e-9 and maxabs(K[3])<1e-9 and err<1e-9
            sixth.append({'L':L,'q':[complex(q).real,complex(q).imag],
                          'f':[complex(f).real,complex(f).imag],
                          'witness':[actual.real,actual.imag],'error':err})
    eigen=[eigen_check(9,0,1j*np.sqrt(3),5),eigen_check(9,1,0,6)]
    out={'seed':20261116,'numpy':np.__version__,'scipy':scipy.__version__,
         'first_order':first,'second_order':seconds,'fourth_order':fourth,
         'fifth_order':fifth,'sixth_order':sixth,'direct_eigenvalue_checks':eigen,
         'parent_identity':parent_checks(rng),'runtime_seconds':time.perf_counter()-start}
    Path(__file__).with_name('checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:(len(v) if isinstance(v,list) else v) for k,v in out.items()},indent=2))


if __name__=='__main__':main()

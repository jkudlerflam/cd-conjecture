"""C19 graded Jacobi-filter cell. Only numpy/scipy and standard library."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import math
from fractions import Fraction
import numpy as np
import scipy.linalg as la

I2=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)

def kron_all(xs):
    a=np.ones((1,1),complex)
    for x in xs: a=np.kron(a,x)
    return a

def modes(n):
    gamma=[]
    for a in range(n):
        gamma += [kron_all([Z]*a+[X]+[I2]*(n-a-1)),
                  kron_all([Z]*a+[Y]+[I2]*(n-a-1))]
    parity=kron_all([Z]*n)
    eta=[1j*g@parity for g in gamma]
    return gamma,eta,parity

def jacobi_integer(j):
    # P_j^(2,1)(2x-1) from its binomial expansion.
    out=[0]*(j+1)
    for l in range(j+1):
        c=math.comb(j+2,l)*math.comb(j+1,j-l)
        for r in range(j-l+1):
            out[l+r]+=c*math.comb(j-l,r)*(-1)**(j-l-r)
    return out

def filters(n):
    coeff=np.zeros((n+1,n),float)
    raw=[]; norms=[]
    for j in range(n):
        q=jacobi_integer(j)
        f=[0]*(n+2)
        for l,v in enumerate(q):
            f[l+1]+=v; f[l+2]-=v
        h=Fraction(j+1,(2*j+4)*(j+3))
        raw.append(f[1:]); norms.append(h)
        coeff[:,j]=np.array(f[1:],float)/math.sqrt(float(h))
    return coeff,raw,norms

def exact_filter_check(n):
    _,raw,hs=filters(n)
    for a in range(n):
        assert sum(raw[a])==0
        for b in range(n):
            value=sum((Fraction(raw[a][e]*raw[b][f],e+f+2)
                       for e in range(n+1) for f in range(n+1)),Fraction(0))
            assert value==(hs[a] if a==b else 0)
    return True

# C26 adaptation of C19. Only the target and an optional visibility multiplier change.
RHO_SCALE=1.0

def primes(count):
    out=[]; x=2
    while len(out)<count:
        if all(x%p for p in out if p*p<=x):out.append(x)
        x+=1
    return out

def sidon_data(n):
    ps=primes(2*n+10)
    q=next(p for p in ps if p>=2*n)
    ids=np.arange(2*n,dtype=int)
    aa=2*q*ids+ids**2%q
    da=int(np.diff(aa).min()); bb=(da+1)//2
    period=2*int(aa[-1])+4*bb+2
    # Preserve the finite resonance margins also at arbitrarily large n.
    rootmax=np.sqrt(ps[2*n-1])
    jitter_scale=min(1e-10,(2*np.pi*da/period)/(8*rootmax))
    if n>=2:
        jitter_scale=min(jitter_scale,(2*np.pi/(3*period))/(16*rootmax))
    jitter=jitter_scale*np.sqrt(np.array(ps[:2*n]))
    points=2*np.pi*(aa+bb+1/6)/period+jitter
    return dict(q=q,A=aa.tolist(),d=da,B=bb,T=period,points=points.tolist(),
                max_jitter=float(jitter.max()),jitter_scale=float(jitter_scale))

def target_parameters(n):
    points=np.array(sidon_data(n)['points'])
    pl,pr=points[:n],points[n:]
    kk=(pr-pl)/2
    theta=-(pr+pl)/2
    t=np.ones(n)
    mu=2*np.cos(kk)
    js=[-0.5*np.array([[np.sin(th),np.cos(th)],[-np.cos(th),np.sin(th)]]) for th in theta]
    return t,theta,mu,np.array([pl,pr]).T,js

def positive_sqrt(a):
    e,u=la.eigh(a)
    assert e.min()>-1e-12
    return (u*np.sqrt(np.maximum(e,0)))@u.conj().T

def make_cell(n):
    N=2**n; k=n+1; r=3
    D=4*math.ceil(max(16,6*k+1)*N/4)
    K=D//4
    gamma,eta,P=modes(n)
    t,theta,mu,momenta,js=target_parameters(n)
    fs,_,_=filters(n)
    p=np.ones(k)/k
    # Syndrome order X,Y,Z. Both A's connect X and Z.
    a1=np.zeros((3,3),complex); a1[0,2]=a1[2,0]=1/math.sqrt(2)
    a2=np.zeros((3,3),complex); a2[0,2]=1j/math.sqrt(2); a2[2,0]=-1j/math.sqrt(2)
    dr=[]; dl=[]
    for e in range(k):
        R=np.zeros((3*N,3*N),complex); L=R.copy()
        for a in range(n):
            for u,A in enumerate([a1,a2]):
                R+=fs[e,a]*np.kron(A,eta[2*a+u])
                op=sum((js[a][u,v]*gamma[2*a+v] for v in range(2)),np.zeros((N,N),complex))
                L-=fs[e,a]*np.kron(A.T,op)
        dr.append(R); dl.append(L)
    largest=max(la.norm(a,2) for a in dr+dl)
    rho=RHO_SCALE/(4*k*largest)
    MR=[p[e]*np.eye(3*N)+rho*dr[e] for e in range(k)]
    ML=[p[e]*np.eye(3*N)+rho*dl[e] for e in range(k)]
    Ms=[la.block_diag(ML[e],MR[e]) for e in range(k)]
    # Physical space: pL,pR,aux. The first 4N auxiliary states hold sL,sR,logical.
    FL=[np.kron(np.kron(f,I2),np.eye(K)) for f in [X,Y,Z]]
    FR=[np.kron(np.kron(I2,f),np.eye(K)) for f in [X,Y,Z]]
    G=np.zeros((D,N),complex)
    for l in range(2):
        for rr in range(2):
            for z in range(N):
                aux=(2*l+rr)*N+z
                G[(2*l+rr)*K+aux,z]=.5
    S=np.concatenate([f@G for f in FL+FR],axis=1)
    extra=[FL[a]@FR[b]@G for a in range(3) for b in range(3)]
    # Four port states for each padding auxiliary state.
    for aux in range(4*N,K):
        for pr in range(4):
            v=np.zeros((D,1),complex);v[pr*K+aux,0]=1;extra.append(v)
    F=np.concatenate([S]+extra,axis=1)
    assert F.shape==(D,D-N)
    auxpar=np.zeros(K)
    logicalpar=np.real(np.diag(P))
    for l in range(2):
        for rr in range(2):
            auxpar[(2*l+rr)*N:(2*l+rr+1)*N]=(-1)**(l+rr)*logicalpar
    for aux in range(4*N,K): auxpar[aux]=(-1)**aux
    parity=np.kron(np.diag(np.kron(Z,Z)).real,auxpar)
    fqpar=np.real(np.einsum('ij,ij->j',F.conj(),parity[:,None]*F))
    sqpar=fqpar[:6*N]
    # W maps the source into energy blocks, and commutes with grading.
    W=np.concatenate([positive_sqrt(m) for m in Ms],axis=0)
    Hq=np.zeros((D-N,D-N),complex)
    es=np.arange(1,k+1,dtype=float)
    for grade in [1,-1]:
        physidx=np.where(np.isclose(fqpar,grade))[0]
        srcidx=np.where(np.isclose(sqpar,grade))[0]
        rows=np.concatenate([e*6*N+srcidx for e in range(k)])
        wg=W[np.ix_(rows,srcidx)]
        d=len(physidx); ns=len(srcidx)
        assert d>=len(rows)
        wp=np.zeros((d,ns),complex);wp[:len(rows)]=wg
        # A complete unitary whose first columns are the specified isometry.
        comp=la.null_space(wp.conj().T)
        U=np.concatenate([wp,comp],axis=1)
        energies=np.r_[np.repeat(es,ns),np.full(d-len(rows),k+1.)]
        h=(U.conj().T*energies)@U
        # F starts with S, so the first physical-grade columns are the source.
        assert np.array_equal(physidx[:ns],srcidx)
        Hq[np.ix_(physidx,physidx)]=h
    H=F@Hq@F.conj().T
    H=(H+H.conj().T)/2
    Kbond=np.zeros((N*N,N*N),complex)
    for a in range(n):
        for u in range(2):
            for v in range(2): Kbond+=js[a][u,v]*np.kron(eta[2*a+u],gamma[2*a+v])
    C=1/(es[:,None]+es[None,:]);scalar=-3*p@C@p
    H2=scalar*np.eye(N*N)+rho*rho*Kbond
    field=sum((-mu[a]/2*kron_all([I2]*a+[Z]+[I2]*(n-a-1)) for a in range(n)),np.zeros((N,N),complex))
    # An abstract even order-lambda^2 on-cell field; annihilates the excited complement.
    field_phys=G@field@G.conj().T
    return dict(n=n,N=N,D=D,K=K,k=k,H=H,G=G,S=S,FL=FL,FR=FR,
                parity=parity,logical_parity=P,ML=ML,MR=MR,Ms=Ms,rho=rho,
                es=es,p=p,filters=fs,H2=H2,Kbond=Kbond,scalar=scalar,
                field=field,field_phys=field_phys,t=t,theta=theta,mu=mu,momenta=momenta)

def direct_second_order(c):
    H,G=c['H'],c['G'];N=c['N'];D=c['D']
    ev,U=la.eigh(H)
    excited=ev>0.5
    ev=ev[excited];U=U[:,excited]
    L=[U.conj().T@f@G for f in c['FL']]
    R=[U.conj().T@f@G for f in c['FR']]
    B=sum((np.kron(R[a],L[a]) for a in range(3)),np.zeros((len(ev)**2,N*N),complex))
    den=(ev[:,None]+ev[None,:]).ravel()
    return -B.conj().T@(B/den[:,None])

def reflection_cell(n, exchange_fields=False):
    """Fully specified uniform completion, no QR or unspecified isometry.

    flag=0 contains the Bell code and source; flags 1,...,k contain W.
    This is a finite matrix cell, not a spatial NN Hamiltonian.
    """
    c=make_cell(n)
    N,k=c['N'],c['k'];K=(k+1)*4*N;D=4*K
    if exchange_fields:
        # d=(1,-1,0,...), p^T C d=1/[2(k+2)].
        g=1/(2*(k+2));p=1/k
        oldrho=c['rho']
        newrho=min(oldrho,np.sqrt(g/(4*k*la.norm(c['field'],2))))
        c['rho']=newrho
        c['ML']=[p*np.eye(3*N)+(newrho/oldrho)*(m-p*np.eye(3*N)) for m in c['ML']]
        c['MR']=[p*np.eye(3*N)+(newrho/oldrho)*(m-p*np.eye(3*N)) for m in c['MR']]
        for e in range(k):
            de=1 if e==0 else (-1 if e==1 else 0)
            c['ML'][e][N:2*N,N:2*N]-=newrho**2/g*de*c['field']
        c['Ms']=[la.block_diag(l,r) for l,r in zip(c['ML'],c['MR'])]
        c['H2']=c['scalar']*np.eye(N*N)+newrho**2*(c['Kbond']+np.kron(np.eye(N),c['field']))
        c['exchange_fields']=True
    FL=[np.kron(np.kron(f,I2),np.eye(K)) for f in [X,Y,Z]]
    FR=[np.kron(np.kron(I2,f),np.eye(K)) for f in [X,Y,Z]]
    def code_at(flag):
        g=np.zeros((D,N),complex)
        for l in range(2):
            for r in range(2):
                for z in range(N):
                    g[(2*l+r)*K+flag*4*N+(2*l+r)*N+z,z]=.5
        return g
    G=code_at(0);S=np.concatenate([f@G for f in FL+FR],axis=1)
    W=np.zeros_like(S)
    for e,m in enumerate(c['Ms']):
        ge=code_at(e+1)
        se=np.concatenate([f@ge for f in FL+FR],axis=1)
        # Use the analytic three-eigenvalue square root on each side.
        root=np.zeros_like(m);p=1/k
        for sl in [slice(0,3*N),slice(3*N,6*N)]:
            dm=(m[sl,sl]-p*np.eye(3*N))/c['rho']
            dm[N:2*N,N:2*N]=0  # Y density field is a separate diagonal block.
            # Explicit eigenvalue radius from row coefficients and t's.
            if sl.start==0:radius=np.sqrt(np.sum(c['filters'][e]**2*c['t']**2)/2)
            else:radius=np.sqrt(2*np.sum(c['filters'][e]**2))
            plus=np.sqrt(p+c['rho']*radius);minus=np.sqrt(p-c['rho']*radius)
            a=(plus-minus)/(2*radius)
            b=(plus+minus-2*np.sqrt(p))/(2*radius**2)
            root[sl,sl]=np.sqrt(p)*np.eye(3*N)+a*dm+b*dm@dm
            offset=sl.start
            yi=slice(offset+N,offset+2*N)
            root[yi,yi]=np.diag(np.sqrt(np.real(np.diag(m[yi,yi]))))
        W+=se@root
    R=np.eye(D)-S@S.conj().T-W@W.conj().T+S@W.conj().T+W@S.conj().T
    energies=np.tile(np.repeat(np.r_[k+1,np.arange(1,k+1)],4*N),4)
    Hdiag=np.diag(energies)-(k+1)*G@G.conj().T
    H=R@Hdiag@R.conj().T
    auxpar=np.tile(np.kron(np.diag(np.kron(Z,Z)).real,np.diag(c['logical_parity']).real),k+1)
    parity=np.kron(np.diag(np.kron(Z,Z)).real,auxpar)
    c.update(D=D,K=K,H=(H+H.conj().T)/2,G=G,S=S,FL=FL,FR=FR,parity=parity,
             field_phys=G@c['field']@G.conj().T,reflection_error=float(la.norm(R@R-np.eye(D))))
    if exchange_fields:c['field_phys']=np.zeros_like(H)
    return c

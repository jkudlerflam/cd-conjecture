import numpy as np, itertools, sys, time
from math import comb, pi
from scipy.sparse.linalg import LinearOperator, eigsh
from scipy.integrate import quad

def harm(n): return sum(1.0/l for l in range(1,n+1))

def bond(N,k):
    # symmetric rank-k rep of SU(N) via occupation numbers
    states=[s for s in itertools.product(range(k+1),repeat=N) if sum(s)==k]
    idx={s:i for i,s in enumerate(states)}; d=len(states)
    E=[[np.zeros((d,d)) for _ in range(N)] for _ in range(N)]
    for s in states:
        for a in range(N):
            for b in range(N):
                if s[b]==0: continue
                t=list(s); amp=np.sqrt(t[b]); t[b]-=1; amp*=np.sqrt(t[a]+1); t[a]+=1
                E[a][b][idx[tuple(t)],idx[s]]+=amp
    X=sum(np.kron(E[a][b],E[b][a]) for a in range(N) for b in range(N))
    w,V=np.linalg.eigh(X)
    # group distinct eigenvalues
    vals=[]; 
    for x in w:
        if not vals or abs(x-vals[-1][0])>1e-8: vals.append([x,1])
        else: vals[-1][1]+=1
    assert len(vals)==k+1, vals
    h=np.zeros((d*d,d*d))
    order=sorted(vals,key=lambda t:-t[0])   # j=0 most symmetric -> largest X
    for j,(x,mult) in enumerate(order):
        sel=np.abs(w-x)<1e-8
        P=V[:,sel]@V[:,sel].T
        h+=2*(harm(k-j)-harm(k))*P
    return d,h,[ (j,2*(harm(k-j)-harm(k)),m) for j,(x,m) in enumerate(order)]

def ground(d,h,L):
    shape=(d,)*L
    def mv(v):
        psi=v.reshape(shape); out=np.zeros_like(psi)
        for i in range(L):
            j=(i+1)%L
            if L==2 and i==1: break
            ax=[i,j]
            t=np.moveaxis(psi,ax,[0,1]); sh=t.shape
            r=(h@t.reshape(d*d,-1)).reshape(sh)
            out+=np.moveaxis(r,[0,1],ax)
        return out.reshape(-1)
    n=d**L
    op=LinearOperator((n,n),matvec=mv,dtype=float)
    rng=np.random.default_rng(1)
    w=eigsh(op,k=1,which='SA',tol=1e-11,v0=rng.standard_normal(n))[0]
    return w[0]

def e0_pred(N,k):
    f=lambda w:(1-np.exp(-k*w))/(2*np.sinh(w/2))*np.sinh((N-1)*w/2)/np.sinh(N*w/2)
    I=quad(f,1e-12,200,limit=500)[0]
    return -2*I

if __name__=="__main__":
    N=int(sys.argv[1]); k=int(sys.argv[2]); Ls=[int(x) for x in sys.argv[3:]]
    d,h,info=bond(N,k)
    print("N,k,d",N,k,d,"bond eigenvalues (j,E,mult):",info)
    c=k*(N*N-1)/(N+k); e0=e0_pred(N,k)
    print("c=",c,"e0 pred=",e0, " v pred=",2*pi/N)
    for L in Ls:
        t=time.time(); E=ground(d,h,L)
        A=-(E-L*e0)*6*L/pi   # = c v effective
        print(f"L={L} E0={E:.10f} E0-L e0={E-L*e0:.8f}  (c v)_eff={A:.5f}  v_eff={A/c:.5f}  v_pred={2*pi/N:.5f} t={time.time()-t:.1f}s",flush=True)

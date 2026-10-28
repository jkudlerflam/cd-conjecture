"""Independent local algebra, Yang-Baxter, dilogarithm and velocity checks.
Run: OPENBLAS_NUM_THREADS=1 python3 fusion_checks.py
Only numpy/scipy plus Python standard library. Writes fusion_checks.json.
"""
import json, math
from itertools import combinations_with_replacement
import numpy as np
from scipy.special import spence

def generators(N,k):
    states=[]
    for inds in combinations_with_replacement(range(N),k):
        states.append(tuple(np.bincount(inds,minlength=N)))
    lookup={x:i for i,x in enumerate(states)}; d=len(states)
    E=np.zeros((N,N,d,d))
    for b in range(N):
      for col,x in enumerate(states):
        if not x[b]: continue
        z=np.array(x); z[b]-=1
        for a in range(N):
          y=z.copy(); y[a]+=1
          E[a,b,lookup[tuple(y)],col]=np.sqrt(x[b]*(z[a]+1))
    return E

def bond(N,k):
    E=generators(N,k); d=E.shape[-1]
    X=sum(np.kron(E[a,b],E[b,a]) for a in range(N) for b in range(N))
    vals,vec=np.linalg.eigh(X)
    eigx=np.array([(k-j)**2-j for j in range(k+1)])
    eh=np.array([-2*sum(1/r for r in range(k-j+1,k+1)) for j in range(k+1)])
    labels=np.argmin(abs(vals[:,None]-eigx[None,:]),axis=1)
    h=(vec*eh[labels])@vec.T
    return X,h,vals,vec,labels,eh

def rogers(x):
    if x>=1-1e-14: return np.pi**2/6
    if x<=1e-14: return 0.
    return spence(1-x)+.5*np.log(x)*np.log1p(-x)

def dilog_c(N,k):
    Q=np.ones((N+1,k+1))
    for r in range(1,N):
      for n in range(1,k):
        Q[r,n]=np.exp(sum(np.log(np.sin(np.pi*(n+j-i)/(N+k))/np.sin(np.pi*(j-i)/(N+k)))
                 for i in range(1,r+1) for j in range(r+1,N+1)))
    total=(N-1)*np.pi**2/6
    for r in range(1,N):
      for n in range(1,k):
        Y=Q[r,n-1]*Q[r,n+1]/(Q[r-1,n]*Q[r+1,n])
        total+=rogers(1/(1+Y))
    return 6*total/np.pi**2

def main():
    out={'local':[],'yang_baxter':[],'dilog':[],'family':[],'velocity':[]}
    for N,k in [(2,1),(2,2),(2,3),(3,2),(7,2),(8,2)]:
      X,h,vals,vec,labels,eh=bond(N,k); d=int(round(np.sqrt(len(h))))
      xpred=np.array([(k-j)**2-j for j in range(k+1)])
      row=dict(N=N,k=k,d=d,X_spectrum_error=float(np.max(np.min(abs(vals[:,None]-xpred),axis=1))),
        h_levels=eh.tolist(), multiplicities=np.bincount(labels,minlength=k+1).tolist(),
        norm=float(np.ptp(eh)/2),trace_per_state=float(np.trace(h)/d**2))
      if k==2: row['polynomial_error']=float(np.max(abs(h-(-X@X/8+3*X/4-np.eye(d*d)))))
      out['local'].append(row)
      if N<=3 and k<=2:
        def R(u):
          r=np.array([np.prod([(l-u)/(l+u) for l in range(k-j+1,k+1)]) for j in range(k+1)])
          return (vec*r[labels])@vec.T
        u,v=.17,.31; I=np.eye(d)
        left=np.kron(R(u),I)@np.kron(I,R(u+v))@np.kron(R(v),I)
        right=np.kron(I,R(v))@np.kron(R(u+v),I)@np.kron(I,R(u))
        out['yang_baxter'].append(dict(N=N,k=k,max_error=float(np.max(abs(left-right))),
            derivative_error=float(np.max(abs((R(1e-6)-R(-1e-6))/2e-6-h)))))
    for N,k in [(2,1),(2,2),(3,2),(7,2),(8,2),(5,5),(8,8),(16,16)]:
      c=k*(N*N-1)/(N+k); M=sum(1/j for j in range(1,k+1))
      out['dilog'].append(dict(N=N,k=k,c_dilog=dilog_c(N,k),c_formula=c))
      out['family'].append(dict(N=N,k=k,d=math.comb(N+k-1,k),c=c,v=2*np.pi/N,norm=M,
             A=2*np.pi*c/N,ratio=c/(N*M)))
    for N in [2,3,8]:
      for r in [1,N-1]:
        lam=3*N; x=2*np.pi*lam/N; a=np.pi*r/N
        # |rho'/rho| directly, no finite-difference loss.
        velocity=2*np.pi/N*np.sinh(x)/(np.cosh(x)-np.cos(a))
        out['velocity'].append(dict(N=N,r=r,rapidity=lam,estimate=velocity,target=2*np.pi/N))
    print(json.dumps(out,indent=2)); open('fusion_checks.json','w').write(json.dumps(out,indent=2))
if __name__=='__main__': main()

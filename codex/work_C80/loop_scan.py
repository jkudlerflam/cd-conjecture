# C80 local copy of C75 source, used as a reference or scalar certificate.
# C75 reuse of C33 scalar-certificate implementation. See PROOF.md.
"""Leading-order loop bands, with tree phases removed before solving.
The real one-cell quadratic term vanishes by particle-hole symmetry and is
set to zero analytically; the raw cancellation residual is retained.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from scipy.linalg import eig
from scipy.optimize import brentq
import itertools,json,time
from structured import Leading,tree,visibility

def margins(roots):
    p=np.sort(roots);twopi=2*np.pi
    def dist(x):return abs((x+np.pi)%twopi-np.pi)
    d2=min(dist(x-y) for i,x in enumerate(p) for y in p[i+1:])
    pairs=list(itertools.combinations_with_replacement(range(len(p)),2))
    vals=sorted(((p[i]+p[j])%twopi,(i,j)) for i,j in pairs)
    diffs=[((vals[(i+1)%len(vals)][0]-vals[i][0])%twopi,vals[i][1],vals[(i+1)%len(vals)][1]) for i in range(len(vals))]
    best=min(diffs)
    return float(d2),float(best[0]),best[1:]

def analyze(n,D,a,eta):
    F,lb=tree(n,D,a);beta=np.exp(lb);T0=-2*visibility(n)*F
    S0=-(T0+T0.T)/2;sc=abs(T0).max();ev=np.linalg.eigvalsh(S0/sc)
    if ev[0]<=0:return dict(n=n,D=D,a=a,eta=eta,status='tree_indefinite')
    td=eta*beta
    T,M,herm=Leading(n,D,a,td=td,loopflux=np.pi).calculate()
    rawreal=np.max(abs(M.real))/sc
    M=1j*M.imag
    M=(M+M.conj().T)/2
    eps=.3*ev[0]*sc
    # Congruence by the definite tree hopping improves the zero polynomial.
    C=np.linalg.cholesky(S0/sc)
    inv=np.linalg.inv(C)
    A=inv@(T/sc)@inv.T;B=inv@((M+eps*np.eye(n))/sc)@inv.T
    I=np.eye(n);Z=np.zeros((n,n))
    zs=eig(np.block([[-B,-A.conj().T],[I,Z]]),np.block([[A,Z],[Z,I]]),right=False)
    unitdev=np.max(abs(abs(zs)-1));roots=np.sort(np.angle(zs))
    def hs(k):return (M+eps*np.eye(n)+T*np.exp(1j*k)+T.conj().T*np.exp(-1j*k))/sc
    # Track roots via eigenvectors, then Newton correction of their Rayleigh quotient.
    residual=[];vel=[];corrections=[]
    for k in roots:
        H=hs(k);w,U=np.linalg.eigh(H);j=np.argmin(abs(w));u=U[:,j]
        v=np.vdot(u,(1j*T*np.exp(1j*k)-1j*T.conj().T*np.exp(-1j*k))@u).real/sc
        residual.append(float(abs(w[j])));vel.append(float(abs(v)))
        corrections.append(float(abs(w[j]/v)))
    d2,d4,witness=margins(roots)
    # Pair each right mover with the left mover of the same ordered band.
    labels=[]
    for k in roots:
        w=np.linalg.eigvalsh(hs(k));labels.append(int(np.argmin(abs(w))))
    centres=[]
    for b in range(n):
        rr=roots[np.array(labels)==b]
        if len(rr)==2:centres.append(float(((sum(rr)/2+np.pi/2)%np.pi)-np.pi/2))
    return dict(n=n,D=D,a=a,eta=eta,td=td,status='numerical',
        T_scale=float(abs(T).max()),tree_herm_min_over_scale=float(ev[0]),
        loop_herm_min_over_scale=float(np.linalg.eigvalsh(-(T+T.conj().T)/2/sc)[0]),
        eps_over_scale=float(eps/sc),mu_over_scale=float(np.linalg.norm(M,2)/sc),
        raw_real_mu_cancellation_residual=float(rawreal),unit_circle_error=float(unitdev),
        max_newton_correction=max(corrections),min_speed_over_scale=min(vel),
        roots=roots.tolist(),centres=centres,delta2=d2,delta4=d4,witness=witness)

if __name__=='__main__':
    out=[]
    for n in range(2,7):
        for eta in [1e-3,1e-5,1e-7]:
            t=time.monotonic();r=analyze(n,100.,2.2,eta);r['seconds']=time.monotonic()-t
            out.append(r);print(json.dumps(r),flush=True)
            json.dump(out,open('loop_results.json','w'),indent=2)

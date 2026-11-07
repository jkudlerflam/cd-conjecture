"""Positivity margin of the pairing argument on the circle (PROOF Section 4).
For each lambda and each out-arc length ell = pi-2Q, compute
  Lam(ell) = J(0) + sup_{x in in} int_out J(x-y) phit(y) dy,  phit=(I-J_out)^{-1}psi, psi=J(.-Q)+J(.+Q),
  Loc(ell) = lower bound on L r0(Q): max(4 r0(Q), (1-omega) r0(Q)/(2 int_Q^{pi/2} r0)),
and report max over ell of Lam/Loc (needs < 1). Floating point, grid N."""
import numpy as np, json, sys
from circle import Circ
def run(lam, N=2048, nell=60):
    M=Circ(np.cosh(lam))
    X=np.linspace(-np.pi/2,np.pi/2,N,endpoint=False)+np.pi/(2*N)
    h=np.pi/N
    Jx=M.J(X)  # J on grid (even)
    def Jd(u):  # J at circular differences via interpolation on fine table
        u=np.mod(u+np.pi/2,np.pi)-np.pi/2
        return np.interp(np.abs(u),Tg,Tv)
    Tg=np.linspace(0,np.pi/2,4001); Tv=M.J(Tg)
    J0=Tv[0]
    worst=(0,None)
    r0=M.r0(X)
    for ell in np.linspace(1e-3,np.pi-1e-3,nell):
        Q=(np.pi-ell)/2
        out=np.abs(X)>Q; Xo=X[out]; Xi=X[~out]
        if out.sum()<2: continue
        Jo=Jd(Xo[:,None]-Xo[None,:])*h
        psi=Jd(Xo-Q)+Jd(Xo+Q)
        phit=np.linalg.solve(np.eye(len(Xo))-Jo,psi)
        Lam=J0+np.max((Jd(Xi[:,None]-Xo[None,:])*h)@phit) if len(Xi) else J0
        omega=np.max(Jo.sum(1))
        r0Q=float(M.r0(np.array([Q]))[0])
        tail=float(np.sum(r0[X>Q])*h)
        Loc=max(4*r0Q,(1-omega)*r0Q/(2*tail))
        if Lam/Loc>worst[0]: worst=(Lam/Loc,float(ell),Lam,Loc,omega)
    return worst
res={}
for lam in [float(a) for a in sys.argv[1].split(',')]:
    w=run(lam); res[lam]=w; print(lam, w, flush=True)

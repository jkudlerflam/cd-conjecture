"""C149: test the uniform barrier Psi = (c1 eta/L) Z_Q + (ce/L) for trial y^+ = quantiles of Z_0 - Psi
(and y^- with -Psi). Reports min_j F_j(y^+)/(2 pi L) * L/eta and max_j F_j(y^-) ... Floating point only."""
import numpy as np, sys
from bethe_ferro import Model, residual, solve
from continuum import Kernel
from trial import continuum_trial

def build(eta, L, n, N=700):
    ker=Kernel(eta); R,P=continuum_trial(eta,L,n,N=N); Q=R.Q; nu=n/L
    # extension of Z_0 beyond Q with s_Q = a1 + int T(u-v) r(v) dv
    yy=np.concatenate([-R.y[::-1],R.y[1:]]); rr=np.concatenate([R.r[::-1],R.r[1:]])
    xs=np.linspace(-1,1,4001)
    def sQ(u):
        # quadrature over [-Q,Q] with breakpoints, simple dense Gauss
        g,w=np.polynomial.legendre.leggauss(40)
        edges=np.unique(np.concatenate([yy, np.clip(u+ eta*np.array([-8,-4,-2,-1,-.5,0,.5,1,2,4,8]),-Q,Q)]))
        tot=0.
        for a,b in zip(edges[:-1],edges[1:]):
            v=(a+b)/2+(b-a)/2*g; tot+=np.sum((b-a)/2*w*ker.T(u-v)*np.interp(np.abs(v),R.y,R.r))
        return ker.a1(u)+tot
    ext=np.concatenate([Q+np.geomspace(1e-6,0.5,60), Q+np.linspace(0.5,25,200)])
    sv=np.array([sQ(u) for u in ext])
    ue=np.concatenate([[Q],ext]); se=np.concatenate([[R.r[-1]],sv])
    cum=np.concatenate([[0],np.cumsum(np.diff(ue)*(se[1:]+se[:-1])/2)])
    def Z0(u):
        a=np.abs(u); s=np.sign(u)
        inside=a<=Q
        out=np.where(inside, P.Z(np.minimum(a,Q)), nu/2+np.interp(a,ue,cum))
        return s*out
    Zfun=lambda u: np.interp(np.abs(u),R.y,R.Z,right=R.Z[-1])
    return R,P,Z0,Zfun

def nodes(Zt, L, n, lo, hi):
    I=np.arange(1,n+1)-(n+1)/2; t=I/L
    a=np.full(n,lo); b=np.full(n,hi)
    for _ in range(200):
        c=(a+b)/2; f=Zt(c)-t; a=np.where(f<0,c,a); b=np.where(f<0,b,c)
    return (a+b)/2

if __name__=='__main__':
    eps=float(sys.argv[1]); L=int(sys.argv[2]); eta=2*np.arcsin(np.sqrt(eps/2)); m=Model(eta)
    c1s=[0.0,0.1,0.2,0.3]; ces=[0.0,0.05,0.1,0.2]
    for n in [L//8,L//4,3*L//8,L//2-8,L//2-2]:
        R,P,Z0,Zf=build(eta,L,n)
        y0=nodes(Z0,L,n,-R.Q-30,R.Q+30); r0=residual(m,y0,L)/(2*np.pi*L)*L/eta
        print(f"n={n} Q={R.Q:.3f} lamQ={eta*L*R.r[-1]:.3g} room*L={eta*(L-2*n)/(2*np.pi):.3f} ZQ(Q)={R.Z[-1]:.2f}  plain: min {r0.min():.3f} max {r0.max():.3f}")
        for c1 in c1s:
            row=[]
            for ce in ces:
                Psi=lambda u: c1*eta/L*Zf(u)+ce/L
                yp=nodes(lambda u: Z0(u)-Psi(u),L,n,-R.Q-30,R.Q+30)
                ym=nodes(lambda u: Z0(u)+Psi(u),L,n,-R.Q-30,R.Q+30)
                Fp=residual(m,yp,L)/(2*np.pi*L)*L/eta; Fm=residual(m,ym,L)/(2*np.pi*L)*L/eta
                row.append(f"{Fp.min():+.3f}/{Fm.max():+.3f}")
            print(f"  c1={c1}: "+"  ".join(f"ce={ce}:{r}" for ce,r in zip(ces,row)),flush=True)
